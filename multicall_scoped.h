#ifndef SCOPED_MULTICALL_H
#define SCOPED_MULTICALL_H

/* multicall_scoped.h
 *
 * 2026, Dana Jacobsen.
 *
 * Based on the idea by Justin Case (justincase / yopmail.com)
 * in https://rt.cpan.org/Public/Bug/Display.html?id=127605.
 *
 * Determines whether we should open a new scope before doing the multicall.
 * It does this by examining the sub's optree.  This gives us high performance
 * in the common case while giving correct behavior in other cases.
 * We consider answering "needs scope" to be the safe thing to do.  Hence we
 * will add a scope if we just aren't sure or the function is too complicated.
 * Keep in mind that a function that complicated is likely to take long enough
 * that our scope won't matter in performance.  We're trying to remove the
 * unnecessary scope for the majority of callbacks which are very simple.
 * SC_MULTICALL is for void context.  SC_MULTICALL_SCALAR stores one return
 * value in the given SV pointer; consume it before the next multicall.
 * SC_MULTICALL_ARRAY appends owned copies of all return values to the given AV.
 *
 * Adds SC_dMULTICALL, SC_PUSH_MULTICALL, the three call variants above, and
 * SC_POP_MULTICALL.
 */

#ifdef dMULTICALL

#ifdef DEBUGGING
#  define SC_ASSERT_GIMME(want) assert(gimme == (want))
#else
#  define SC_ASSERT_GIMME(want) ((void)0)
#endif

#define SC_dMULTICALL              \
  dMULTICALL;                      \
  bool multicall_needs_scope_ = 0; \
  SV *multicall_scalar_ = NULL;

#define SC_PUSH_MULTICALL(cv)      \
  STMT_START {                     \
    PUSH_MULTICALL(cv);            \
    multicall_needs_scope_ = sc_multicall_cv_needs_scope(aTHX_ (cv)); \
    if (multicall_needs_scope_ && gimme == G_SCALAR) { \
      multicall_scalar_ = newSV(0); \
      SAVEFREESV(multicall_scalar_); \
    }                              \
  } STMT_END

#define SC_MULTICALL               \
  STMT_START {                     \
    SC_ASSERT_GIMME(G_VOID);       \
    if (multicall_needs_scope_)    \
      ENTER;                       \
    MULTICALL;                     \
    if (multicall_needs_scope_)    \
      LEAVE;                       \
  } STMT_END

/* LEAVE may release a lexical SV returned by the callback. */
#define SC_MULTICALL_SCALAR(sv)    \
  STMT_START {                     \
    SC_ASSERT_GIMME(G_SCALAR);     \
    if (multicall_needs_scope_) {  \
      ENTER;                       \
      MULTICALL;                   \
      sv_setsv(multicall_scalar_, *PL_stack_sp); \
      LEAVE;                       \
      (sv) = multicall_scalar_;    \
    } else {                       \
      MULTICALL;                   \
      (sv) = *PL_stack_sp;         \
    }                              \
  } STMT_END

/* Copy list results before LEAVE can release lexical return values. */
#define SC_MULTICALL_ARRAY(av)     \
  STMT_START {                     \
    SSize_t sc_before_ix_ = (SSize_t)(PL_stack_sp - PL_stack_base); \
    SSize_t sc_i_, sc_nret_;       \
    AV *sc_av_ = (av);             \
    SC_ASSERT_GIMME(G_ARRAY);      \
    if (multicall_needs_scope_) {  \
      /* pp_leavesub expects stack slot zero to remain the undef sentinel. */ \
      if (sc_before_ix_ <= 0)      \
        *PL_stack_base = &PL_sv_undef; \
      ENTER;                       \
    }                              \
    MULTICALL;                     \
    sc_nret_ = (SSize_t)(PL_stack_sp - PL_stack_base) - sc_before_ix_; \
    for (sc_i_ = 1; sc_i_ <= sc_nret_; sc_i_++) \
      av_push(sc_av_, newSVsv(PL_stack_base[sc_before_ix_ + sc_i_])); \
    if (multicall_needs_scope_)    \
      LEAVE;                       \
    PL_stack_sp = PL_stack_base + sc_before_ix_; \
  } STMT_END

#define SC_POP_MULTICALL  POP_MULTICALL


/* Roll our own tiny PERL_VERSION_GE */
#ifndef PERL_SUBVERSION
#  ifdef SUBVERSION
#    define PERL_SUBVERSION SUBVERSION
#  else
#    define PERL_SUBVERSION 0
#  endif
#endif
#define SC_PERL_AT_LEAST(v, s) \
  ((PERL_VERSION > (v)) || (PERL_VERSION == (v) && PERL_SUBVERSION >= (s)))


/* Returns 0 if we see no reason to wrap this sub inside it's own scope.
   Returns 1 if we need to because of locals created.
   Returns 1 if it's too complicated (long, infinite loop, deep branches) */

#define SC_MAX_BRANCHES 8
#define SC_MAX_NOPS     500
static bool sc_multicall_cv_needs_scope(pTHX_ const CV *cv) {
  OP *o = CvSTART(cv);
  size_t nops = 0;
  OP *branches[SC_MAX_BRANCHES];
  int nbranch = 0;
  for (;; o = o->op_next) {
    if (!o) {
      if (nbranch > 0) { o = branches[--nbranch]; continue; }
      break;
    }
    /* printf("   %s\n",PL_op_name[o->op_type]); */
    if (++nops > SC_MAX_NOPS) return 1; /* Too big */
    switch (o->op_type) {
      case OP_PADSV:  case OP_PADAV:  case OP_PADHV:
      case OP_ANONCODE:
#if SC_PERL_AT_LEAST(17,6)
      case OP_PADRANGE:
#endif
#if SC_PERL_AT_LEAST(27,6)
      case OP_MULTICONCAT:  /* This could hide a PADSV -- we don't know */
#endif
#if SC_PERL_AT_LEAST(37,3)
      case OP_PADSV_STORE:
#endif
#if SC_PERL_AT_LEAST(43,3)
      case OP_MULTIPARAM:
#endif
        return 1;

      case OP_AND:  case OP_OR:  case OP_COND_EXPR:
      case OP_ANDASSIGN:  case OP_ORASSIGN:
#if SC_PERL_AT_LEAST(9,0)
      case OP_DOR:
      case OP_DORASSIGN:
#endif
        if (nbranch >= SC_MAX_BRANCHES) return 1; /* Too deep */
        branches[nbranch++] = cLOGOPx(o)->op_other;
        break;
      case OP_LEAVESUB:
        if (nbranch > 0) { o = branches[--nbranch]; continue; }
        break;
    }
  }
  return 0;
}

#endif

#endif
