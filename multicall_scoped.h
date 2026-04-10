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
 *
 * Adds SC_ versions of dMULTICALL, PUSH_MULTICALL, and MULTICALL.
 */

#ifdef dMULTICALL

#define SC_dMULTICALL              \
  dMULTICALL;                      \
  bool multicall_needs_scope_ = 0;

#define SC_PUSH_MULTICALL(cv)      \
  STMT_START {                     \
    PUSH_MULTICALL(cv);            \
    multicall_needs_scope_ = sc_multicall_cv_needs_scope(aTHX_ (cv)); \
  } STMT_END

#define SC_MULTICALL               \
  STMT_START {                     \
    if (multicall_needs_scope_)    \
      ENTER;                       \
    MULTICALL;                     \
    if (multicall_needs_scope_)    \
      LEAVE;                       \
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
