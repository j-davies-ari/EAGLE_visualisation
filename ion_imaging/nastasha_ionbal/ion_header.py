import eagle_constants_and_units as c


column = {'elec'    : 0,
          'HI'      : 1,
          'HII'     : 2,
          'Hm'      : 3,
          'HeI'     : 4,
          'HeII'    : 5,
          'HeIII'   : 6,
          'CI'      : 7,
          'CII'     : 8,
          'CIII'    : 9,
          'CIV'     : 10,
          'CV'      : 11,
          'CVI'     : 12,
          'CVII'    : 13,
          'Cm'      : 14,
          'NI'      : 15,
          'NII'     : 16,
          'NIII'    : 17,
          'NIV'     : 18,
          'NV'      : 19,
          'NVI'     : 20,
          'NVII'    : 21,
          'NVIII'   : 22,
          'OI'      : 23,
          'OII'     : 24,
          'OIII'    : 25,
          'OIV'     : 26,
          'OV'      : 27,
          'OVI'     : 28,
          'OVII'    : 29,
          'OVIII'   : 30,
          'OIX'     : 31,
          'Om'      : 32,
          'NeI'     : 33,
          'NeII'    : 34,
          'NeIII'   : 35,
          'NeIV'    : 36,
          'NeV'     : 37,
          'NeVI'    : 38,
          'NeVII'   : 39,
          'NeVIII'  : 40,
          'NeIX'    : 41,
          'NeX'     : 42,
          'NeXI'    : 43,
          'MgI'     : 44,
          'MgII'    : 45,
          'MgIII'   : 46,
          'MgIV'    : 47,
          'MgV'     : 48,
          'MgVI'    : 49,
          'MgVII'   : 50,
          'MgVIII'  : 51,
          'MgIX'    : 52,
          'MgX'     : 53,
          'MgXI'    : 54,
          'MgXII'   : 55,
          'MgXIII'  : 56,
          'SiI'     : 57,
          'SiII'    : 58,
          'SiIII'   : 59,
          'SiIV'    : 60,
          'SiV'     : 61,
          'SiVI'    : 62,
          'SiVII'   : 63,
          'SiVIII'  : 64,
          'SiIX'    : 65,
          'SiX'     : 66,
          'SiXI'    : 67,
          'SiXII'   : 68,
          'SiXIII'  : 69,
          'SiXIV'   : 70,
          'SiXV'    : 71,
          'SI'      : 72,
          'SII'     : 73,
          'SIII'    : 74,
          'SIV'     : 75,
          'SV'      : 76,
          'SVI'     : 77,
          'SVII'    : 78,
          'SVIII'   : 79,
          'SIX'     : 80,
          'SX'      : 81,
          'SXI'     : 82,
          'SXII'    : 83,
          'SXIII'   : 84,
          'SXIV'    : 85,
          'SXV'     : 86,
          'SXVI'    : 87,
          'SXVII'   : 88,
          'CaI'     : 89,
          'CaII'    : 90,
          'CaIII'   : 91,
          'CaIV'    : 92,
          'CaV'     : 93,
          'CaVI'    : 94,
          'CaVII'   : 95,
          'CaVIII'  : 96,
          'CaIX'    : 97,
          'CaX'     : 98,
          'CaXI'    : 99,
          'CaXII'   : 100,
          'CaXIII'  : 101,
          'CaXIV'   : 102,
          'CaXV'    : 103,
          'CaXVI'   : 104,
          'CaXVII'  : 105,
          'CaXVIII' : 106,
          'CaXIX'   : 107,
          'CaXX'    : 108,
          'CaXXI'   : 109,
          'FeI'     : 110,
          'FeII'    : 111,
          'FeIII'   : 112,
          'FeIV'    : 113,
          'FeV'     : 114,
          'FeVI'    : 115,
          'FeVII'   : 116,
          'FeVIII'  : 117,
          'FeIX'    : 118,
          'FeX'     : 119,
          'FeXI'    : 120,
          'FeXII'   : 121,
          'FeXIII'  : 122,
          'FeXIV'   : 123,
          'FeXV'    : 124,
          'FeXVI'   : 125,
          'FeXVII'  : 126,
          'FeXVIII' : 127,
          'FeXIX'   : 128,
          'FeXX'    : 129,
          'FeXXI'   : 130,
          'FeXXII'  : 131,
          'FeXXIII' : 132,
          'FeXXIV'  : 133,
          'FeXXV'   : 134,
          'FeXXVI'  : 135,
          'FeXXVII' : 136,
          'H2'      : 137,
          'H2p'     : 138,
          'H3p'     : 139,
          'OH'      : 140,
          'H2O'     : 141,
          'C2'      : 142,
          'O2'      : 143,
          'HCOp'    : 144,
          'CH'      : 145,
          'CH2'     : 146,
          'CH3p'    : 147,
          'CO'      : 148,
          'CHp'     : 149,
          'CH2p'    : 150,
          'OHp'     : 151,
          'H2Op'    : 152,
          'H3Op'    : 153,
          'COp'     : 154,
          'HOCp'    : 155,
          'O2p'     : 156}

atomw = {'Hydrogen' : c.atomw_H,
         'HI'       : c.atomw_H,
         'HII'      : c.atomw_H,
         'Helium'   : c.atomw_He,
         'HeI'      : c.atomw_He,
         'HeII'     : c.atomw_He,
         'HeIII'    : c.atomw_He,
         'Carbon'   : c.atomw_C,
         'CI'       : c.atomw_C,
         'CII'      : c.atomw_C,
         'CIII'     : c.atomw_C,
         'CIV'      : c.atomw_C,
         'CV'       : c.atomw_C,
         'CVI'      : c.atomw_C,
         'CVII'     : c.atomw_C,
	 'Nitrogen' : c.atomw_N,
         'NI'       : c.atomw_N,
         'NII'      : c.atomw_N,
         'NIII'     : c.atomw_N,
         'NIV'      : c.atomw_N,
         'NV'       : c.atomw_N,
         'NVI'      : c.atomw_N,
         'NVII'     : c.atomw_N,
         'NVIII'    : c.atomw_N,
         'Oxygen'   : c.atomw_O,
         'OI'       : c.atomw_O,
         'OII'      : c.atomw_O,
         'OIII'     : c.atomw_O,
         'OIV'      : c.atomw_O,
         'OV'       : c.atomw_O,
         'OVI'      : c.atomw_O,
         'OVII'     : c.atomw_O,
         'OVIII'    : c.atomw_O,
         'OIX'      : c.atomw_O,
         'Neon'     : c.atomw_Ne,
         'NeI'      : c.atomw_Ne,
         'NeII'     : c.atomw_Ne,
         'NeIII'    : c.atomw_Ne,
         'NeIV'     : c.atomw_Ne,
         'NeV'      : c.atomw_Ne,
         'NeVI'     : c.atomw_Ne,
         'NeVII'    : c.atomw_Ne,
         'NeVIII'   : c.atomw_Ne,
         'NeIX'     : c.atomw_Ne,
         'NeX'      : c.atomw_Ne,
         'NeXI'     : c.atomw_Ne,
         'Magnesium': c.atomw_Mg,
         'MgI'      : c.atomw_Mg,
         'MgII'     : c.atomw_Mg,
         'MgIII'    : c.atomw_Mg,
         'MgIV'     : c.atomw_Mg,
         'MgV'      : c.atomw_Mg,
         'MgVI'     : c.atomw_Mg,
         'MgVII'    : c.atomw_Mg,
         'MgVIII'   : c.atomw_Mg,
         'MgIX'     : c.atomw_Mg,
         'MgX'      : c.atomw_Mg,
         'MgXI'     : c.atomw_Mg,
         'MgXII'    : c.atomw_Mg,
         'MgXIII'   : c.atomw_Mg,
         'Silicon'  : c.atomw_Si,
         'SiI'      : c.atomw_Si,
         'SiII'     : c.atomw_Si,
         'SiIII'    : c.atomw_Si,
         'SiIV'     : c.atomw_Si,
         'SiV'      : c.atomw_Si,
         'SiVI'     : c.atomw_Si,
         'SiVII'    : c.atomw_Si,
         'SiVIII'   : c.atomw_Si,
         'SiIX'     : c.atomw_Si,
         'SiX'      : c.atomw_Si,
         'SiXI'     : c.atomw_Si,
         'SiXII'    : c.atomw_Si,
         'SiXIII'   : c.atomw_Si,
         'SiXIV'    : c.atomw_Si,
         'SiXV'     : c.atomw_Si,
         'Sulfur'   : c.atomw_S,
         'SI'       : c.atomw_S,
         'SII'      : c.atomw_S,
         'SIII'     : c.atomw_S,
         'SIV'      : c.atomw_S,
         'SV'       : c.atomw_S,
         'SVI'      : c.atomw_S,
         'SVII'     : c.atomw_S,
         'SVIII'    : c.atomw_S,
         'SIX'      : c.atomw_S,
         'SX'       : c.atomw_S,
         'SXI'      : c.atomw_S,
         'SXII'     : c.atomw_S,
         'SXIII'    : c.atomw_S,
         'SXIV'     : c.atomw_S,
         'SXV'      : c.atomw_S,
         'SXVI'     : c.atomw_S,
         'SXVII'    : c.atomw_S,
         'Calcium'  : c.atomw_Ca,
         'CaI'      : c.atomw_Ca,
         'CaII'     : c.atomw_Ca,
         'CaIII'    : c.atomw_Ca,
         'CaIV'     : c.atomw_Ca,
         'CaV'      : c.atomw_Ca,
         'CaVI'     : c.atomw_Ca,
         'CaVII'    : c.atomw_Ca,
         'CaVIII'   : c.atomw_Ca,
         'CaIX'     : c.atomw_Ca,
         'CaX'      : c.atomw_Ca,
         'CaXI'     : c.atomw_Ca,
         'CaXII'    : c.atomw_Ca,
         'CaXIII'   : c.atomw_Ca,
         'CaXIV'    : c.atomw_Ca,
         'CaXV'     : c.atomw_Ca,
         'CaXVI'    : c.atomw_Ca,
         'CaXVII'   : c.atomw_Ca,
         'CaXVIII'  : c.atomw_Ca,
         'CaXIX'    : c.atomw_Ca,
         'CaXX'     : c.atomw_Ca,
         'CaXXI'    : c.atomw_Ca,
         'Iron'     : c.atomw_Fe,
         'FeI'      : c.atomw_Fe,
         'FeII'     : c.atomw_Fe,
         'FeIII'    : c.atomw_Fe,
         'FeIV'     : c.atomw_Fe,
         'FeV'      : c.atomw_Fe,
         'FeVI'     : c.atomw_Fe,
         'FeVII'    : c.atomw_Fe,
         'FeVIII'   : c.atomw_Fe,
         'FeIX'     : c.atomw_Fe,
         'FeX'      : c.atomw_Fe,
         'FeXI'     : c.atomw_Fe,
         'FeXII'    : c.atomw_Fe,
         'FeXIII'   : c.atomw_Fe,
         'FeXIV'    : c.atomw_Fe,
         'FeXV'     : c.atomw_Fe,
         'FeXVI'    : c.atomw_Fe,
         'FeXVII'   : c.atomw_Fe,
         'FeXVIII'  : c.atomw_Fe,
         'FeXIX'    : c.atomw_Fe,
         'FeXX'     : c.atomw_Fe,
         'FeXXI'    : c.atomw_Fe,
         'FeXXII'   : c.atomw_Fe,
         'FeXXIII'  : c.atomw_Fe,
         'FeXXIV'   : c.atomw_Fe,
         'FeXXV'    : c.atomw_Fe,
         'FeXXVI'   : c.atomw_Fe,
         'FeXXVII'  : c.atomw_Fe}
