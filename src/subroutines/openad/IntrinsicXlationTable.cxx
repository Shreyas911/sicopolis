// -*-Mode: C++;-*-
// $Header: /Volumes/cvsrep/developer/OpenADFortTk/src/lib/support/IntrinsicXlationTable.cxx,v 1.28 2006/05/12 16:12:22 utke Exp $

#include <stdlib.h>  
#include <algorithm> 


#include "IntrinsicXlationTable.h"
#include "Open64IRInterface/wn_attr.h"
#include "Diagnostics.h"

namespace fortTkSupport { 

  IntrinsicXlationTable::Entry IntrinsicXlationTable::ourTable[] = {

    // -------------------------------------------------------
    // WHIRL calls, expressions and intrinsics that correspond to XAIF
    // Intrinsics.  
    //
    // For OPR_INTRINSIC_OP, the WHIRL string is the INTRINSIC_basename()
    // (cf. wintrinsic.h, wutil.cxx)
    // -------------------------------------------------------

    // Common mathematical functions
    { { WNExpr,     OPR_NEG,          NULL,      1, 0 }, { XAIFIntrin, "minus_scal",                NULL,            1 } },
    { { WNExpr,     OPR_ADD,          NULL,      2, 0 }, { XAIFIntrin, "add_scal_scal",             NULL,            2 } },
    { { WNExpr,     OPR_SUB,          NULL,      2, 0 }, { XAIFIntrin, "sub_scal_scal",             NULL,            2 } },
    { { WNExpr,     OPR_MPY,          NULL,      2, 0 }, { XAIFIntrin, "mul_scal_scal",             NULL,            2 } },
    { { WNExpr,     OPR_DIV,          NULL,      2, 0 }, { XAIFIntrin, "div_scal_scal",             NULL,            2 } }, 
    { { WNCall,     OPR_CALL,         "SQRT",    1, 0 }, { XAIFIntrin, "sqrt_scal",                 "0_SQRT",        1 } },
    { { WNCall,     OPR_CALL,         "DSQRT",   1, 0 }, { XAIFIntrin, "sqrt_scal",                 "1_DSQRT",       1 } },
    { { WNExpr,     OPR_SQRT,         NULL,      1, 0 }, { XAIFIntrin, "sqrt_scal",                 "2_SQRT",        1 } }, 
    { { WNCall,     OPR_CALL,         "SUM",     1, 0 }, { XAIFIntrin, "sum_arr",                   NULL,            1 } }, 
    // modulo/remainder
    { { WNExpr,     OPR_MOD,          NULL,      2, 0 }, { XAIFIntrin, "bogus_modulo_scal_scal",    "0_MODULO",      2 } }, 
    { { WNCall,     OPR_CALL,         "MODULO",  2, 0 }, { XAIFIntrin, "bogus_modulo_scal_scal",    "1_MODULO",      2 } }, 
    { { WNExpr,     OPR_REM,          NULL,      2, 0 }, { XAIFIntrin, "mod_scal_scal",             "0_MOD",         2 } },
    { { WNIntrinOp, OPR_INTRINSIC_OP, "AMOD",    2, 0 }, { XAIFIntrin, "mod_scal_scal",             "1_AMOD",        2 } },
    { { WNIntrinOp, OPR_INTRINSIC_OP, "DMOD",    2, 0 }, { XAIFIntrin, "mod_scal_scal",             "2_DMOD",        2 } }, 
    { { WNCall,     OPR_CALL,         "MOD",     2, 0 }, { XAIFIntrin, "mod_scal_scal",             "3_MOD",         2 } }, 
    // trigonometric
    { { WNCall,     OPR_CALL,         "SIN",     1, 0 }, { XAIFIntrin, "sin_scal",                  "0_SIN",         1 } },
    { { WNCall,     OPR_CALL,         "DSIN",    1, 0 }, { XAIFIntrin, "sin_scal",                  "1_DSIN",        1 } },
    { { WNCall,     OPR_CALL,         "COS",     1, 0 }, { XAIFIntrin, "cos_scal",                  "0_COS",         1 } },
    { { WNCall,     OPR_CALL,         "DCOS",    1, 0 }, { XAIFIntrin, "cos_scal",                  "1_DCOS",        1 } },
    { { WNCall,     OPR_CALL,         "TAN",     1, 0 }, { XAIFIntrin, "tan_scal",                  "0_TAN",         1 } },
    { { WNCall,     OPR_CALL,         "OAD_TAN", 1, 0 }, { XAIFIntrin, "oad_tan_scal",              "0_TAN",         1 } },
    { { WNCall,     OPR_CALL,         "DTAN",    1, 0 }, { XAIFIntrin, "tan_scal",                  "1_DTAN",        1 } }, 
    { { WNCall,     OPR_CALL,         "ASIN",    1, 0 }, { XAIFIntrin, "arcsin_scal",               "0_ASIN",        1 } },
    { { WNCall,     OPR_CALL,         "ACOS",    1, 0 }, { XAIFIntrin, "arccos_scal",               "0_ACOS",        1 } },
    { { WNCall,     OPR_CALL,         "ATAN",    1, 0 }, { XAIFIntrin, "arctan_scal",               "0_ATAN",        1 } }, 
    { { WNCall,     OPR_CALL,         "ATAN2",   2, 0 }, { XAIFIntrin, "arctan_scal_scal",          "0_ATAN2",       2 } }, 
    { { WNCall,     OPR_CALL,         "SINH",    1, 0 }, { XAIFIntrin, "sinh_scal",                 "0_SINH",        1 } },
    { { WNCall,     OPR_CALL,         "DSINH",   1, 0 }, { XAIFIntrin, "sinh_scal",                 "1_DSINH",       1 } },
    { { WNCall,     OPR_CALL,         "COSH",    1, 0 }, { XAIFIntrin, "cosh_scal",                 "0_COSH",        1 } },
    { { WNCall,     OPR_CALL,         "DCOSH",   1, 0 }, { XAIFIntrin, "cosh_scal",                 "1_DCOSH",       1 } },
    { { WNCall,     OPR_CALL,         "TANH",    1, 0 }, { XAIFIntrin, "tanh_scal",                 "0_TANH",        1 } },
    { { WNCall,     OPR_CALL,         "DTANH",   1, 0 }, { XAIFIntrin, "tanh_scal",                 "1_DTANH",       1 } },
    // exp/log
    { { WNCall,     OPR_CALL,         "EXP",     1, 0 }, { XAIFIntrin, "exp_scal",                  "0_EXP",         1 } },
    { { WNCall,     OPR_CALL,         "DEXP",    1, 0 }, { XAIFIntrin, "exp_scal",                  "1_DEXP",        1 } },
    { { WNCall,     OPR_CALL,         "LOG",     1, 0 }, { XAIFIntrin, "ln_scal",                   "0_LOG",         1 } },
    { { WNCall,     OPR_CALL,         "DLOG",    1, 0 }, { XAIFIntrin, "ln_scal",                   "1_DLOG",        1 } }, 
    { { WNCall,     OPR_CALL,         "ALOG",    1, 0 }, { XAIFIntrin, "ln_scal",                   "2_ALOG",        1 } }, 
    { { WNCall,     OPR_CALL,         "LOG10",   1, 0 }, { XAIFIntrin, "log10",                     NULL,            1 } }, 
    { { WNIntrinOp, OPR_INTRINSIC_OP, "EXPEXPR", 2, 0 }, { XAIFIntrin, "pow_scal_scal",             NULL,            2 } },
    // string operations
    { { WNIntrinOp, OPR_INTRINSIC_OP, "CEQEXPR", 2, 0 }, { XAIFIntrin, "string_eq_scal_scal",       NULL,            2 } },
    { { WNIntrinOp, OPR_INTRINSIC_OP, "CNEEXPR", 2, 0 }, { XAIFIntrin, "string_ne_scal_scal",       NULL,            2 } },
    { { WNIntrinOp, OPR_INTRINSIC_OP, "CGEEXPR", 2, 0 }, { XAIFIntrin, "string_ge_scal_scal",       NULL,            2 } },
    { { WNIntrinOp, OPR_INTRINSIC_OP, "CGTEXPR", 2, 0 }, { XAIFIntrin, "string_gt_scal_scal",       NULL,            2 } },
    { { WNIntrinOp, OPR_INTRINSIC_OP, "CLEEXPR", 2, 0 }, { XAIFIntrin, "string_le_scal_scal",       NULL,            2 } },
    { { WNIntrinOp, OPR_INTRINSIC_OP, "CLTEXPR", 2, 0 }, { XAIFIntrin, "string_lt_scal_scal",       NULL,            2 } },
    { { WNIntrinOp, OPR_INTRINSIC_OP, "LEN",     1, 0 }, { XAIFIntrin, "len",                       "0_LEN_OP",      1 } },
    { { WNCall,     OPR_CALL,         "INDEX",   3, 1 }, { XAIFIntrin, "index",                     "0_INDEX_CALL",  3 } },
    { { WNCall,     OPR_CALL,         "LEN",     1, 0 }, { XAIFIntrin, "len",                       "1_LEN_CALL",    1 } },
    { { WNCall,     OPR_CALL,         "LEN_TRIM",1, 0 }, { XAIFIntrin, "len_trim",                  NULL,            1 } },
    { { WNCall,     OPR_CALL,         "TRIM",    1, 0 }, { XAIFIntrin, "trim",                      NULL,            1 } },
    { { WNCall,     OPR_CALL,         "SCAN",    2, 0 }, { XAIFIntrin, "scan",                      NULL,            2 } },
    { { WNCall,     OPR_CALL,         "ICHAR",   1, 0 }, { XAIFIntrin, "ichar",                     NULL,            1 } },
    // Rounding, conversion
    { { WNExpr,     OPR_ABS,          NULL,      1, 0 }, { XAIFIntrin, "abs_scal",                  "0_ABS",         1 } },
    { { WNCall,     OPR_CALL,         "ABS",     1, 0 }, { XAIFIntrin, "abs_scal",                  "1_ABS",         1 } },
    { { WNCall,     OPR_CALL,         "DABS",    1, 0 }, { XAIFIntrin, "abs_scal",                  "2_DABS",        1 } },
    { { WNCall,     OPR_CALL,         "IABS",    1, 0 }, { XAIFIntrin, "iabs_scal",                 "0_IABS",        1 } }, 
    { { WNCall,     OPR_CALL,         "SIGN",    2, 0 }, { XAIFIntrin, "sign_scal_scal",            "0_SIGN",        2 } }, 
    { { WNCall,     OPR_CALL,         "DSIGN",   2, 0 }, { XAIFIntrin, "sign_scal_scal",            "1_SIGN",        2 } }, 
    { { WNExpr,     OPR_RND,          NULL,      1, 0 }, { XAIFIntrin, "bogus_round_scal",          NULL,            1 } },
    { { WNExpr,     OPR_TRUNC,        NULL,      1, 0 }, { XAIFIntrin, "int_scal",                  "0_INT",         1 } },
    { { WNCall,     OPR_CALL,         "INT",     1, 0 }, { XAIFIntrin, "int_scal",                  "1_INT",         1 } },
    { { WNCall,     OPR_CALL,         "NINT",    1, 0 }, { XAIFIntrin, "nint_scal",                 NULL,            1 } },
    { { WNCall,     OPR_CALL,         "TRANSFER",2, 0 }, { XAIFIntrin, "transfer",                  NULL,            2 } },
    { { WNExpr,     OPR_CEIL,         NULL,      1, 0 }, { XAIFIntrin, "bogus_ceil_scal",           NULL,            1 } },
    { { WNExpr,     OPR_FLOOR,        NULL,      1, 0 }, { XAIFIntrin, "bogus_floor_scal",          NULL,            1 } }, 
    { { WNCall,     OPR_CALL,         "REAL",    1, 0 }, { XAIFIntrin, "real_scal",                 "0_REAL",        1 } },
    { { WNCall,     OPR_CALL,         "FLOAT",   1, 0 }, { XAIFIntrin, "real_scal",                 "1_REAL",        1 } },
    { { WNCall,     OPR_CALL,         "DBLE",    1, 0 }, { XAIFIntrin, "real_scal",                 "2_REAL",        1 } }, 
    { { WNCall,     OPR_CALL,         "AIMAG",   1, 0 }, { XAIFIntrin, "imag_scal",                 NULL,            1 } }, 
    { { WNCall,     OPR_CALL,         "TRANSPOSE",1, 0 }, { XAIFIntrin, "transpose_arr",            NULL,            1 } }, 
    { { WNCall,     OPR_CALL,         "RESHAPE", 2, 0 }, { XAIFIntrin, "reshape_arr",               NULL,            2 } }, 
    { { WNExpr,     OPR_COMPLEX,      NULL,      2, 0 }, { XAIFIntrin, "complex_scal",              NULL,            2 } }, 
    // Logical (and bitwise logical) operations 
    { { WNExpr,     OPR_BNOT,         NULL,      1, 0 }, { XAIFBoolOp, "b_not",                     NULL,            1 } },
    { { WNExpr,     OPR_BAND,         NULL,      2, 0 }, { XAIFBoolOp, "b_and",                     NULL,            2 } },
    { { WNExpr,     OPR_BIOR,         NULL,      2, 0 }, { XAIFBoolOp, "b_or",                      NULL,            2 } },
    { { WNExpr,     OPR_BXOR,         NULL,      2, 0 }, { XAIFBoolOp, "b_xor",                     NULL,            2 } }, 
    { { WNExpr,     OPR_LNOT,         NULL,      1, 0 }, { XAIFBoolOp, "not",                       NULL,            1 } },
    { { WNExpr,     OPR_LAND,         NULL,      2, 0 }, { XAIFBoolOp, "and",                       NULL,            2 } },
    { { WNExpr,     OPR_LIOR,         NULL,      2, 0 }, { XAIFBoolOp, "or",                        NULL,            2 } },
    { { WNExpr,     OPR_CAND,         NULL,      2, 0 }, { XAIFBoolOp, "bogus_cand_scal_scal",      NULL,            2 } },
    { { WNExpr,     OPR_CIOR,         NULL,      2, 0 }, { XAIFBoolOp, "bogus_cor_scal_scal",       NULL,            2 } }, 
    { { WNExpr,     OPR_EQ,           NULL,      2, 0 }, { XAIFBoolOp, "equal",                     NULL,            2 } },
    { { WNExpr,     OPR_NE,           NULL,      2, 0 }, { XAIFBoolOp, "not_equal",                 NULL,            2 } },
    { { WNExpr,     OPR_GT,           NULL,      2, 0 }, { XAIFBoolOp, "greater_than",              NULL,            2 } },
    { { WNExpr,     OPR_GE,           NULL,      2, 0 }, { XAIFBoolOp, "greater_or_equal",          NULL,            2 } },
    { { WNExpr,     OPR_LT,           NULL,      2, 0 }, { XAIFBoolOp, "less_than",                 NULL,            2 } },
    { { WNExpr,     OPR_LE,           NULL,      2, 0 }, { XAIFBoolOp, "less_or_equal",             NULL,            2 } }, 
    { { WNCall,     OPR_CALL,         "ANY",     1, 0 }, { XAIFIntrin, "any",                       NULL,            1 } },
    // Misc.
    { { WNCall,     OPR_CALL,         "MAXLOC",  1, 0 }, { XAIFIntrin, "maxloc_arr",                NULL,            1 } },
    { { WNCall,     OPR_CALL,         "MINLOC",  1, 0 }, { XAIFIntrin, "minloc_arr",                NULL,            1 } },
    { { WNCall,     OPR_CALL,         "LBOUND",  2, 0 }, { XAIFIntrin, "lbound",                    NULL,            2 } },
    { { WNCall,     OPR_CALL,         "UBOUND",  2, 0 }, { XAIFIntrin, "ubound",                    NULL,            2 } },
    { { WNCall,     OPR_CALL,         "SIZE",    2, 1 }, { XAIFIntrin, "size",                      NULL,            2 } },
    { { WNCall,     OPR_CALL,         "SHAPE",   1, 1 }, { XAIFIntrin, "shape_arr",                 NULL,            1 } },
    { { WNIntrinOp, OPR_INTRINSIC_OP, "F90INDEX",2, 0 }, { XAIFIntrin, "index",                     "1_INDEX_OP",    2 } },
    { { WNExpr,     OPR_SHL,           NULL,     2, 0 }, { XAIFIntrin, "bogus_shl_scal_scal",       NULL,            2 } },
    { { WNExpr,     OPR_ASHR,          NULL,     2, 0 }, { XAIFIntrin, "bogus_ashr_scal_scal",      NULL,            2 } }, 
    { { WNCall,     OPR_CALL,         "PRESENT", 1, 0 }, { XAIFIntrin, "present",                   NULL,            1 } },
    { { WNCall,     OPR_CALL,         "ASSOCIATED", 1, 0 }, { XAIFIntrin, "associated",             NULL,            1 } },
    { { WNCall,     OPR_CALL,         "ALLOCATED", 1, 0 }, { XAIFIntrin, "allocated",               NULL,            1 } },
    { { WNCall,     OPR_CALL,         "_ALLOCATE", 1, 0 }, { XAIFIntrin, "allocate",                NULL,            1 } },
    { { WNIntrinOp, OPR_NULLIFY,      NULL     , 1, 0 }, { XAIFIntrin, "nullify",                   NULL,            1 } },
    { { WNCall,     OPR_CALL,         "_DEALLOCATE", 1, 0 }, { XAIFIntrin, "deallocate",            NULL,            1 } },
    // max/min etc. are turned into special subroutine calls by the canonicalizer except for integer expressions;
    // nc stands for not-canonicalized 
    { { WNExpr,     OPR_MAX,           NULL,     2, 0 }, { XAIFIntrin, "nc_max_scal_scal",          NULL,            2 } },
    { { WNCall,     OPR_CALL,         "MAXVAL",  1, 0 }, { XAIFIntrin, "nc_maxval",                 NULL,            1 } },
    { { WNExpr,     OPR_MIN,           NULL,     2, 0 }, { XAIFIntrin, "nc_min_scal_scal",          NULL,            2 } },
    { { WNCall,     OPR_CALL,         "MINVAL",  1, 0 }, { XAIFIntrin, "nc_minval",                 NULL,            1 } }
  
  };

  unsigned int IntrinsicXlationTable::ourTableSize = 
    (sizeof(IntrinsicXlationTable::ourTable) / sizeof(IntrinsicXlationTable::Entry));

  const std::string IntrinsicXlationTable::toString(const WNOprClass& oprcl) {
    std::string retStr;
    switch (oprcl) {
    case WNCall:     
      retStr="WNCall";
      break;
    case WNIntrinCall:      
      retStr="WNIntrinCall";
      break;
    case WNIntrinOp:      
      retStr="WNIntrinOp";
      break;
    case WNExpr:     
      retStr="WNExpr";
      break;
    default:
      FORTTK_DIE("IntrinsicXlationTable::toString: unknown WNOprClass " << (int)oprcl);
      break;
    }
    return retStr; // should never reach
  }

  const std::string IntrinsicXlationTable::toString(const XAIFOpr& opr) {
    std::string retStr;
    switch (opr) {
    case XAIFSubCall:        
      retStr="XAIFSubCall";
      break;
    case XAIFFuncCall:          
      retStr="XAIFFuncCall";
      break;
    case XAIFIntrin:            
      retStr="XAIFIntrinsic";
      break;
    case XAIFBoolOp:            
      retStr="XAIFBoolOp";
      break;
    default:
      FORTTK_DIE("IntrinsicXlationTable::toString: unknown XAIFOpr "<< (int)opr);
      break;
    }
    return retStr; 
  }

  void IntrinsicXlationTable::WHIRLInfo::dump(std::ostream& os) const {
    os << "{ " 
       << toString(oprcl).c_str() << ", " 
       << OPERATOR_name(opr)      << ", "
       << ((name) ? name : "<null>") << ", "
       << numop 
       << " }";
  }

  void IntrinsicXlationTable::WHIRLInfo::ddump() const {
    dump(std::cerr);
  }

  void IntrinsicXlationTable::XAIFInfo::dump(std::ostream& os) const {
    os << "{ " 
       << toString(opr).c_str() << ", " 
       << ((name) ? name : "<null>") << ", " 
       << ((key) ? key : "<null>") << ", " 
       << numop 
       << " }";
  }

  void IntrinsicXlationTable::XAIFInfo::ddump() const {
    dump(std::cerr);
  }

  IntrinsicXlationTable::LtSortedTable::LtSortedTable(TableType aTableType, bool ignoreXaifKey) : 
    myTableType(aTableType), 
    myIgnoreXaifKeyFlag(ignoreXaifKey) { 
  }

  IntrinsicXlationTable::LtSortedTable::~LtSortedTable() { 
  }

  bool IntrinsicXlationTable::LtSortedTable:: operator()(const Entry* e1, 
							 const Entry* e2) const {
    if (myTableType == W2X) {
      // 1. whirl_info.opr is the primary sorting key
      if (e1->whirl_info.opr == e2->whirl_info.opr) {
      
	// Either 1) both whirl_info.oprcl will be equal or 2) one will
	// be 'WNOprClass_UNKNOWN' (the search item)
	WNOprClass cl = (e1->whirl_info.oprcl == WNOprClass_UNKNOWN) ? 
	  e2->whirl_info.oprcl : e1->whirl_info.oprcl;

	switch (cl) {
	case WNCall:
	case WNIntrinCall:
	case WNIntrinOp:
	  // 2. whirl_info.name is the secondary sorting key
	  return (strcmp(e1->whirl_info.name, e2->whirl_info.name) < 0);
	case WNExpr:
	  // 2. There is no secondary sorting key
	  return false; // e1 and e2 are equal
	default:
	  FORTTK_DIE("Internal IntrinsicXlationTable error: Unknown WNOprClass: "
		     << cl);
	}
      } 
      else {
	return (e1->whirl_info.opr < e2->whirl_info.opr);
      }
    } 
    else if (myTableType == X2W) {
      // 1. xaif_info.opr is the primary sorting key
      if (e1->xaif_info.opr == e2->xaif_info.opr) {
	// 2. xaif_info.name is the secondary sorting key
	int cmp = (strcmp(e1->xaif_info.name, e2->xaif_info.name));
	if (!myIgnoreXaifKeyFlag && cmp == 0) {
	  // 3. xaif_info.key, if available, is the tertiary sorting key
	  if (e1->xaif_info.key && e2->xaif_info.key) {
	    return (strcmp(e1->xaif_info.key, e2->xaif_info.key) < 0);
	  }
	  else if (e1->xaif_info.key) /* && !e2->xaif_info.key */ {
	    return false; // e1 > e2
	  }
	  else if (e2->xaif_info.key) /* && !e1->xaif_info.key */ {
	    return true; // e1 < e2
	  }
	  // fall-through
	} 
	return (cmp < 0);
      } 
      else {
	return (e1->xaif_info.opr < e2->xaif_info.opr);
      }
    } 
    else {
      FORTTK_DIE("Internal IntrinsicXlationTable error: Unknown TableType: " 
		 << myTableType);
    }
  }

  IntrinsicXlationTable::IntrinsicXlationTable(const TableType& aTableType) : 
    myTableType(aTableType), 
    mySortedTable(SortedTable(ourTableSize)) {
    // Initialize it
    for (unsigned int i = 0; i < ourTableSize; ++i) {
      Entry* e = &ourTable[i];
      mySortedTable[i] = e;
    }
    // Sort it ascendingly
    std::sort(mySortedTable.begin(), 
	      mySortedTable.end(), 
	      LtSortedTable(myTableType));
  }

  IntrinsicXlationTable::~IntrinsicXlationTable() 
  { 
  }

  IntrinsicXlationTable::XAIFInfoPair::XAIFInfoPair(bool aBoolean,
						    const IntrinsicXlationTable::XAIFInfo& anXAIFInfo) :
    first(aBoolean),
    second(anXAIFInfo) { 
  }

  IntrinsicXlationTable::XAIFInfoPair IntrinsicXlationTable::findXAIFInfo(OPERATOR opr, 
									  const char* name,
									  bool mustFind) {
    if (myTableType != W2X) 
      FORTTK_DIE("IntrinsicXlationTable::findXAIFInfo: wrong TableType"); 
    static Entry toFind;
    toFind.whirl_info.oprcl = WNOprClass_UNKNOWN;
    toFind.whirl_info.opr = opr;
    toFind.whirl_info.name = name;
    LtSortedTable lt(myTableType);
    SortedTableIt it = std::lower_bound(mySortedTable.begin(), 
					mySortedTable.end(),
					toFind, lt);
    if (it != mySortedTable.end() && !lt(*it, toFind) && !lt(toFind, *it)) {
      return XAIFInfoPair(true,const_cast<const XAIFInfo&>((*it)->xaif_info));
    } 
    if (mustFind)
      FORTTK_DIE("IntrinsicXlationTable::findXAIFInfo: unknown opr="
		 << OPERATOR_name(opr)
		 << ",name="
		 << name);
    return XAIFInfoPair(false,toFind.xaif_info);
  }

  IntrinsicXlationTable::WHIRLInfo* 
  IntrinsicXlationTable::findWHIRLInfo(XAIFOpr opr, 
				       const char* name, 
				       const char* key) {
    if (myTableType != X2W)
      FORTTK_DIE("IntrinsicXlationTable::findWHIRLInfo: wrong TableType"); 
    static Entry toFind;  
    toFind.xaif_info.opr = opr;
    toFind.xaif_info.name = name;
    toFind.xaif_info.key = (key && key[0] == '\0') ? NULL : key;
    bool ignoreKey = (toFind.xaif_info.key == NULL);
    LtSortedTable lt(myTableType, ignoreKey);
    SortedTableIt it = std::lower_bound(mySortedTable.begin(), 
					mySortedTable.end(), 
					toFind, 
					lt);
    if (it != mySortedTable.end() && !lt(*it, toFind) && !lt(toFind, *it)) {
      Entry* foundItem = *it; // may be multiple matches if key is NULL
      return &(foundItem->whirl_info);
    } 
    FORTTK_DIE("IntrinsicXlationTable::findWHIRLInfo: unknown opr="
	       << toString(opr).c_str()
	       << ",name="
	       << name
	       << ",key="
	       << key);
    return NULL;
  }

  class PrintEntry : public std::unary_function<IntrinsicXlationTable::Entry*, void> {
  public:
    PrintEntry(std::ostream& anOstream, IntrinsicXlationTable::TableType aTableType) : 
      myOstream(anOstream), 
      myTableType(aTableType) { 
    }
    ~PrintEntry() { 
    }
    void operator() (IntrinsicXlationTable::Entry* anEntry) {
      if (myTableType == IntrinsicXlationTable::W2X) {
	anEntry->whirl_info.dump(myOstream); 
	myOstream << " ";
	anEntry->xaif_info.dump(myOstream);  
	myOstream << std::endl;
      } 
      else if (myTableType == IntrinsicXlationTable::X2W) {
	anEntry->xaif_info.dump(myOstream);  
	myOstream << " ";
	anEntry->whirl_info.dump(myOstream); 
	myOstream << std::endl;
      } 
      else {
	FORTTK_DIE("Internal IntrinsicXlationTable error: Unknown TableType: " 
		   << myTableType);
      }
    }
  private:
    std::ostream& myOstream;
    IntrinsicXlationTable::TableType myTableType;
  };

  void IntrinsicXlationTable::dump(std::ostream& os) const {
    os << "Begin Intrinsic Table\n";
    std::for_each(mySortedTable.begin(), 
		  mySortedTable.end(), 
		  PrintEntry(os, 
			     myTableType));
    os << "End Intrinsic Table" << std::endl;
  }

  void IntrinsicXlationTable::ddump() const {
    dump(std::cerr);
  }

}
