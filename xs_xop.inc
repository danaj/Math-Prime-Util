#if MPU_HAS_CUSTOM_OPS

typedef OP* (*xop_ppfunc_t)(pTHX);
typedef struct {
  const char* cv_name;
  const char* xop_name;
  const char* xop_desc;
  xop_ppfunc_t ppfunc;
  U8 nargs;
  XOP xop;
} xop_registration_t;
static xop_registration_t* xop_find_registration(GV *namegv);

static OP* xop_call_checker_exact_arity(pTHX_ OP *entersubop, GV *namegv, SV *ckobj) {
  xop_registration_t *xopreg = xop_find_registration(namegv);
  OP *pushop, *cur, *arg_first, *arg_last, *cvop, *newop;
  U8 i;

  PERL_UNUSED_ARG(ckobj);
  if (!xopreg) return entersubop;

  pushop = cUNOPx(entersubop)->op_first;
  if (!pushop) return entersubop;
  if (!OpHAS_SIBLING(pushop))
    pushop = cUNOPx(pushop)->op_first;
  if (!pushop) return entersubop;
  if (pushop->op_type == OP_NULL && cUNOPx(pushop)->op_first != NULL)
    pushop = cUNOPx(pushop)->op_first;
  if (pushop->op_type != OP_PUSHMARK) return entersubop;

  arg_first = arg_last = NULL;
  cur = OpSIBLING(pushop);
  for (i = 0; i < xopreg->nargs; i++) {
    if (!cur) return entersubop;
    /* Previously we were rejecting (returning entersubop) anything with $_.
     * Later this was limited to 5.18, and now removed entirely as I am no
     * longer finding issues.  It seems like this was a trigger that showed
     * bugs in our stack manipulation later.  Unlike simpler XS code, we have
     * lot of back and forth XS<->Perl<->XS<->... callbacks possible. */
    if (i == 0) arg_first = cur;
    arg_last = cur;
    cur = OpSIBLING(cur);
  }
  cvop = cur;
  if (!cvop) return entersubop;
  if (OpSIBLING(cvop) != NULL) return entersubop;
  if (cvop->op_type != OP_GV && cvop->op_type != OP_RV2CV) {
    if (cvop->op_type == OP_NULL && cUNOPx(cvop)->op_first != NULL) {
      OP *cvkid = cUNOPx(cvop)->op_first;
      if (cvkid->op_type != OP_GV && cvkid->op_type != OP_RV2CV) return entersubop;
    } else {
      return entersubop;
    }
  }

  /* The custom op inherits its caller's context.  Preserve the scalar
   * argument context supplied by the original prototypes, even when the
   * custom op itself is used in void context. */
  cur = arg_first;
  for (i = 0; i < xopreg->nargs; i++) {
    op_contextualize(cur, G_SCALAR);
    cur = OpSIBLING(cur);
  }

  OpMORESIB_set(pushop, cvop);
  if (xopreg->nargs > 0)
    OpLASTSIB_set(arg_last, NULL);

  if (xopreg->nargs == 0)
    newop = newOP(OP_NULL, 0);
  else if (xopreg->nargs == 1)
    newop = newUNOP(OP_NULL, 0, arg_first);
  else if (xopreg->nargs == 2)
    newop = newBINOP(OP_NULL, 0, arg_first, arg_last);
  else
    newop = newLISTOP(OP_NULL, 0, arg_first, arg_last);
  newop->op_type = OP_CUSTOM;
  newop->op_ppaddr = xopreg->ppfunc;

  op_free(entersubop);
  return newop;
}

/* This is for xop GMP calls that return one bigint-like scalar. */
static void xop_canonicalize_result(pTHX) {
  dSP;
  SV* out = xs_to_canonical(aTHX_ TOPs);
  SPAGAIN;
  SETs(out);
  PUTBACK;
}

/* Calling call_sv from a custom pp op can let the nested runloop execute the
 * caller's continuation before call_sv returns.  Hide the current op during
 * the callback; SAVEOP restores it on both return and croak. */
static vcall_used_t xop_call_backend(pTHX_ const char* name, int nargs, int minversion) {
  OP *curop = PL_op;
  vcall_used_t used = VCALL_USED_NONE;
  ENTER;
  SAVEOP();
  PL_op = NULL;
  (void)_vcallsubn(aTHX_ G_SCALAR, VCALL_PP|VCALL_GMP, name, nargs, minversion, &used);
  LEAVE;
  assert(PL_op == curop);
  return used;
}

static OP* xop_dispatch_binary(pTHX_ const char* name, int minversion) {
  dSP;
  vcall_used_t used;
  PUTBACK;
  used = xop_call_backend(aTHX_ name, 2, minversion);
  if (used == VCALL_USED_GMP) xop_canonicalize_result(aTHX);
  SPAGAIN;
  RETURN;
}

static OP* xop_dispatch_unary(pTHX_ const char* name, int minversion) {
  dSP;
  vcall_used_t used;
  PUTBACK;
  used = xop_call_backend(aTHX_ name, 1, minversion);
  if (used == VCALL_USED_GMP) xop_canonicalize_result(aTHX);
  SPAGAIN;
  RETURN;
}

static OP* pp_irand_custom(pTHX) {
  dSP;
  dMY_CXT;
  PUSHs(sv_2mortal(newSVuv(irand32(MY_CXT.randcxt))));
  RETURN;
}

static OP* pp_irand64_custom(pTHX) {
  dSP;
  dMY_CXT;
#if BITS_PER_WORD == 32
  PUSHs(sv_2mortal(newSVuv(irand32(MY_CXT.randcxt))));
#else
  PUSHs(sv_2mortal(newSVuv(irand64(MY_CXT.randcxt))));
#endif
  RETURN;
}

static OP* pp_drand_custom(pTHX) {
  dSP;
  dMY_CXT;
  PUSHs(sv_2mortal(newSVnv(drand64(MY_CXT.randcxt))));
  RETURN;
}

static OP* pp_validate_integer_custom(pTHX) {
  dSP;
  xs_validate_integer_inplace(aTHX_ TOPs, IFLAG_ANY);
  SPAGAIN;
  SETs(&PL_sv_yes);
  RETURN;
}
static OP* pp_validate_integer_nonneg_custom(pTHX) {
  dSP;
  xs_validate_integer_inplace(aTHX_ TOPs, IFLAG_NONNEG);
  SPAGAIN;
  SETs(&PL_sv_yes);
  RETURN;
}
static OP* pp_validate_integer_positive_custom(pTHX) {
  dSP;
  xs_validate_integer_inplace(aTHX_ TOPs, IFLAG_POS);
  SPAGAIN;
  SETs(&PL_sv_yes);
  RETURN;
}
static OP* pp_validate_integer_abs_custom(pTHX) {
  dSP;
  xs_validate_integer_inplace(aTHX_ TOPs, IFLAG_ABS);
  SPAGAIN;
  SETs(&PL_sv_yes);
  RETURN;
}

static OP* pp_addint_custom_common(pTHX_ int opix, const char *opname, int minversion) {
  dSP;
  SV *sva = SP[-1], *svb = SP[0];
  SV *slowret;
  int astatus, bstatus, is_uv;
  UV uvret;
  IV ivret;

  if (addint_try_native_result(aTHX_ opix, sva, svb, opname,
                               &astatus, &bstatus, &is_uv, &uvret, &ivret)) {
    SP -= 1;
    SETs(sv_2mortal(is_uv ? newSVuv(uvret) : newSViv(ivret)));
    RETURN;
  }

  slowret = addint_try_slow_result(aTHX_ opix, sva, svb, astatus, bstatus, opname);
  if (slowret != NULL) {
    SPAGAIN;
    SP -= 1;
    SETs(slowret);
    RETURN;
  }

  return xop_dispatch_binary(aTHX_ opname, minversion);
}

static OP* pp_addint_custom(pTHX)  { return pp_addint_custom_common(aTHX_ 0, "addint", 52); }
static OP* pp_subint_custom(pTHX)  { return pp_addint_custom_common(aTHX_ 1, "subint", 52); }
static OP* pp_mulint_custom(pTHX)  { return pp_addint_custom_common(aTHX_ 2, "mulint", 52); }
static OP* pp_divint_custom(pTHX)  { return pp_addint_custom_common(aTHX_ 3, "divint", 52); }
static OP* pp_modint_custom(pTHX)  { return pp_addint_custom_common(aTHX_ 4, "modint", 52); }
static OP* pp_cdivint_custom(pTHX) { return pp_addint_custom_common(aTHX_ 5, "cdivint", 53); }
static OP* pp_powint_custom(pTHX)  { return pp_addint_custom_common(aTHX_ 7, "powint", 54); }

static OP* pp_add1int_custom_common(pTHX_ int opix, const char *opname) {
  dSP;
  SV *svn = TOPs;
  SV *slowret;
  int status, is_uv;
  UV uvret;
  IV ivret;

  if (add1_try_native_result(aTHX_ opix, svn, &status, &is_uv, &uvret, &ivret)) {
    SETs(sv_2mortal(is_uv ? newSVuv(uvret) : newSViv(ivret)));
    RETURN;
  }

  slowret = add1_try_slow_result(aTHX_ opix, svn);
  if (slowret != NULL) {
    SPAGAIN;
    SETs(slowret);
    RETURN;
  }

  return xop_dispatch_unary(aTHX_ opname, 53);
}

static OP* pp_add1int_custom(pTHX) { return pp_add1int_custom_common(aTHX_ 0, "add1int"); }
static OP* pp_sub1int_custom(pTHX) { return pp_add1int_custom_common(aTHX_ 1, "sub1int"); }

static OP* pp_signint_custom_common(pTHX_ int opix, const char *opname) {
  dSP;
  SV *svret, *svn = TOPs;
  int status, sign, isodd, ret;
  UV n;
  const char* s;
  STRLEN len;

  status = _validate_and_set(&n, aTHX_ svn, IFLAG_ANY);
  if (status == 0) {
    s = SvPV_nomg(svn, len);
    if (len == 0 || s == 0) croak("%s: invalid non-empty input", opname);
    sign = (s[0] == '-') ? -1 : (s[0] == '0') ? 0 : 1;
    isodd = (s[len-1] == '1' || s[len-1] == '3' || s[len-1] == '5' || s[len-1] == '7' || s[len-1] == '9');
  } else {
    sign = (status == -1) ? -1 : (n == 0) ? 0 : 1;
    isodd = n & 1;
  }
  ret = (opix==0) ? sign : (opix==1) ? isodd : !isodd;
  svret = _sv_const_int(aTHX_ ret);
  SETs(svret);
  RETURN;
}
static OP* pp_signint_custom(pTHX) { return pp_signint_custom_common(aTHX_ 0, "signint"); }
static OP* pp_is_odd_custom(pTHX)  { return pp_signint_custom_common(aTHX_ 1, "is_odd"); }
static OP* pp_is_even_custom(pTHX) { return pp_signint_custom_common(aTHX_ 2, "is_even"); }

static OP* pp_is_square_custom(pTHX) {
  dSP;
  UV n;
  int ret;
  int status = _validate_and_set(&n, aTHX_ TOPs, IFLAG_ANY);

  if (status != 0) {
    SETs(_sv_const_int(aTHX_ status == 1 && is_perfect_square(n)));
    RETURN;
  }
  if (xs_sv_is_perfect_square(aTHX_ TOPs, &ret)) {
    SETs(_sv_const_int(aTHX_ ret));
    RETURN;
  }
  return xop_dispatch_unary(aTHX_ "is_square", 47);
}

static OP* pp_cmpint_custom(pTHX) {
  dSP;
  SV *sva = SP[-1], *svb = SP[0];
  int astatus, bstatus, ret = 0;
  UV a, b;

  astatus = _validate_and_set(&a, aTHX_ sva, IFLAG_ANY);
  bstatus = _validate_and_set(&b, aTHX_ svb, IFLAG_ANY);
  if (astatus != 0 && bstatus != 0) {
    if      (astatus > bstatus) ret = 1;
    else if (astatus < bstatus) ret = -1;
    else if (a == b)            ret = 0;
    else                        ret = ((astatus == 1 && a > b) || (astatus == -1 && (IV)a > (IV)b)) ? 1 : -1;
  } else {
    STRLEN alen, blen;
    char *aptr, *bptr;
    aptr = SvPV(sva, alen);
    bptr = SvPV(svb, blen);
    ret = strint_cmp(aptr, alen, bptr, blen);
  }
  SP -= 1;
  SETs(_sv_const_int(aTHX_ ret));
  RETURN;
}

static OP* pp_kronecker_custom(pTHX) {
  dSP;
  SV *sva = SP[-1], *svb = SP[0];
  int astatus, bstatus, ret;
  UV a, b;

  astatus = _validate_and_set(&a, aTHX_ sva, IFLAG_ANY);
  bstatus = _validate_and_set(&b, aTHX_ svb, IFLAG_ANY);
  if (astatus != 0 && bstatus != 0) {
    if (bstatus == 1)
      ret = (astatus == 1) ? kronecker_uu(a, b) : kronecker_su((IV)a, b);
    else
      ret = (astatus == 1) ? kronecker_uu(a, neg_iv(b)) : -kronecker_su((IV)a, neg_iv(b));
    SP -= 1;
    SETs(_sv_const_int(aTHX_ ret));
    RETURN;
  }
  return xop_dispatch_binary(aTHX_ "kronecker", 17);
}

static xop_registration_t xop_registrations[] = {
  { "Math::Prime::Util::_validate_integer", "_validate_integer", "_validate_integer custom op", pp_validate_integer_custom, 1, {0} },
  { "Math::Prime::Util::_validate_integer_nonneg", "_validate_integer_nonneg", "_validate_integer_nonneg custom op", pp_validate_integer_nonneg_custom, 1, {0} },
  { "Math::Prime::Util::_validate_integer_positive", "_validate_integer_positive", "_validate_integer_positive custom op", pp_validate_integer_positive_custom, 1, {0} },
  { "Math::Prime::Util::_validate_integer_abs", "_validate_integer_abs", "_validate_integer_abs custom op", pp_validate_integer_abs_custom, 1, {0} },
  { "Math::Prime::Util::irand", "irand", "irand custom op", pp_irand_custom, 0, {0} },
  { "Math::Prime::Util::irand64", "irand64", "irand64 custom op", pp_irand64_custom, 0, {0} },
  { "Math::Prime::Util::drand", "drand", "drand custom op", pp_drand_custom, 0, {0} },
  { "Math::Prime::Util::addint", "addint", "addint custom op", pp_addint_custom, 2, {0} },
  { "Math::Prime::Util::subint", "subint", "subint custom op", pp_subint_custom, 2, {0} },
  { "Math::Prime::Util::add1int", "add1int", "add1int custom op", pp_add1int_custom, 1, {0} },
  { "Math::Prime::Util::sub1int", "sub1int", "sub1int custom op", pp_sub1int_custom, 1, {0} },
  { "Math::Prime::Util::mulint", "mulint", "mulint custom op", pp_mulint_custom, 2, {0} },
  { "Math::Prime::Util::divint", "divint", "divint custom op", pp_divint_custom, 2, {0} },
  { "Math::Prime::Util::modint", "modint", "modint custom op", pp_modint_custom, 2, {0} },
  { "Math::Prime::Util::cdivint", "cdivint", "cdivint custom op", pp_cdivint_custom, 2, {0} },
  { "Math::Prime::Util::powint", "powint", "powint custom op", pp_powint_custom, 2, {0} },
  { "Math::Prime::Util::is_odd", "is_odd", "is_odd custom op", pp_is_odd_custom, 1, {0} },
  { "Math::Prime::Util::is_even", "is_even", "is_even custom op", pp_is_even_custom, 1, {0} },
  { "Math::Prime::Util::is_square", "is_square", "is_square custom op", pp_is_square_custom, 1, {0} },
  { "Math::Prime::Util::cmpint", "cmpint", "cmpint custom op", pp_cmpint_custom, 2, {0} },
  { "Math::Prime::Util::kronecker", "kronecker", "kronecker custom op", pp_kronecker_custom, 2, {0} },
  { "Math::Prime::Util::signint", "signint", "signint custom op", pp_signint_custom, 1, {0} },

  /* We considered absint, negint, sqrtint.  Returning a possible bigint from
   * those operations is more difficult than it looks. */
};

static xop_registration_t* xop_find_registration(GV *namegv) {
  const char *name;
  U32 i;
  if (!namegv || !isGV(namegv)) return NULL;
  name = GvNAME(namegv);
  if (!name) return NULL;
  for (i = 0; i < (U32)(sizeof(xop_registrations)/sizeof(xop_registrations[0])); i++) {
    if (strEQ(name, xop_registrations[i].xop_name))
      return &xop_registrations[i];
  }
  return NULL;
}

#endif
