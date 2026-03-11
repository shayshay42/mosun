# Manual explicit ODE/rule core for the mosun model.
# Runtime code includes this file directly.

const DIRECT_STATE_COUNT = 131
const DIRECT_PARAM_COUNT = 107
const DIRECT_VALUE_COUNT = 238

const DIRECT_STATE_NAMES = [
    "actTtiss",
    "Btiss",
    "TDBc_ugperkg",
    "actTpb",
    "drugBtisskill",
    "drugTtissact",
    "restTtiss",
    "TDBp_ugperkg",
    "restTpb",
    "Bpb",
    "TDBc_ugperml",
    "TDBt_ugperml",
    "drugTpbact",
    "drugBpbkill",
    "restTtiss_perml",
    "actTtiss_perml",
    "Btiss_perml",
    "restTpb_perml",
    "actTpb_perml",
    "Bpb_perml",
    "BTrratio_pb",
    "BTrratio_tiss",
    "TaBratio_pb",
    "TaBratio_tiss",
    "Tafraction_pb",
    "totTpb_perml",
    "totTtiss_perml",
    "act0Ttiss_perml",
    "act0Tpb_perml",
    "Tafraction_tiss",
    "BAFF",
    "Baffconsumption",
    "act0Tpb",
    "act0Ttiss",
    "injection_effect",
    "Btiss2",
    "act0Tpb_1",
    "Bpb_1",
    "restTpb_1",
    "actTpb_1",
    "act0Ttiss2",
    "restTtiss2",
    "TDBt2_ugperml",
    "drugTtissact2",
    "drugBtisskill2",
    "actTtiss2",
    "Tafraction_tiss2",
    "act0Ttiss2_perml",
    "totTtiss2_perml",
    "TaBratio_tiss2",
    "BTrratio_tiss2",
    "Btiss2_perml",
    "actTtiss2_perml",
    "restTtiss2_perml",
    "RTXc_ugperkg",
    "RTXp_ugperkg",
    "RTXc_ugperml",
    "RTXt_ugperml",
    "RTXt2_ugperml",
    "drug_effect",
    "B1920tiss3",
    "B19no20tiss3",
    "restTtiss3",
    "act0Ttiss3",
    "actTtiss3",
    "drugTtissact3",
    "TDBt3_ugperml",
    "B19Trratio_tiss3",
    "B19tiss3",
    "restTtiss3_perml",
    "act0Ttiss3_perml",
    "actTtiss3_perml",
    "B19tiss3_perml",
    "Tafraction_tiss3",
    "totTtiss3_perml",
    "TaB19ratio_tiss3",
    "B19TtotRatio_tiss3",
    "RTXt3_ugperml",
    "drugBtisskill3",
    "TaB1920ratio_tiss3",
    "B1920tiss3_perml",
    "B1920Trratio_tiss3",
    "B19no20tiss3_perml",
    "Blinc_ug",
    "Blinc_ngperml",
    "BlinTpbact",
    "BlinBpbkill",
    "BlinTtissact",
    "BlinBtisskill",
    "Blint_ngperml",
    "Blint2_ngperml",
    "BlinBtisskill2",
    "BlinTtissact2",
    "BlinBtisskill3",
    "Blint3_ngperml",
    "BlinTtissact3",
    "Bpb_norm",
    "Tafraction_pb_init",
    "totTtiss",
    "totTtiss2",
    "totTtiss3",
    "restTtumor",
    "actTtumor",
    "Btumor",
    "act0Ttumor",
    "drugTtumoract",
    "drugBtumorkill",
    "BlinTtumoract",
    "BlinBtumorkill",
    "TDBtumor_ugperml",
    "restTtumor_perml",
    "act0Ttumor_perml",
    "actTtumor_perml",
    "Btumor_perml",
    "Tafraction_tumor",
    "BTrratio_tumor",
    "TaBratio_tumor",
    "totTtumor",
    "totTtumor_perml",
    "Blintumor_ngperml",
    "RTXtumor_ugperml",
    "BTtotRatio_tiss",
    "BTtotRatio_tiss2",
    "TDBsc_ugperkg",
    "TDBc_ugperml_AUC",
    "IL6pb",
    "IL6tiss",
    "IL6tiss2",
    "IL6tiss3",
    "IL6tumor",
    "IL6combo",
]

const DIRECT_PARAM_NAMES = [
    "VmB",
    "KmTB_kill",
    "KdrugB",
    "kBapop",
    "kBprolif",
    "kTprolif",
    "kBkill",
    "KdrugactT",
    "VmT",
    "KmBT_act",
    "kTaexit",
    "kTact",
    "fTadeact",
    "fTap",
    "Cl_tdb",
    "Cld_tdb",
    "Vc_tdb",
    "Vp_tdb",
    "kTgen",
    "kTaapop",
    "KTrp",
    "Trpbo_perml",
    "fTaprolif",
    "fBprolif",
    "Bpbo_perml",
    "fBexit",
    "KBp",
    "Kp",
    "Vpb",
    "Vtissue",
    "fTrapop",
    "Trpbref_perml",
    "nkill",
    "KTrp2",
    "Vm_tdb",
    "Km_tdb",
    "BW",
    "act0on",
    "kIL6prod",
    "fa0",
    "fAICD",
    "kBAFFprod",
    "Bpbref_perml",
    "thBAFF",
    "thalfIL6",
    "fTa0deact",
    "fTa0apop",
    "BAFFo",
    "fBAFFo",
    "tinjhalf",
    "finj",
    "KBp2",
    "Vtissue2",
    "Kp2",
    "tissue2on",
    "Vc_rtx",
    "Vp_rtx",
    "Cl_rtx",
    "Cld_rtx",
    "Vm_rtx",
    "Km_rtx",
    "Kmkill_rtx",
    "fTgenbl",
    "depleteTpb",
    "depleteBpb",
    "kTrexit",
    "PKflag",
    "VPid",
    "end_time",
    "fdrug",
    "Vtissue3",
    "KTrp3",
    "tissue3on",
    "Kp3",
    "kBtiss3exit",
    "KBp3",
    "B19no20_B1920_ratio",
    "fvalidation",
    "fapop_v24",
    "fBtissue3_v1",
    "kBmat_kBapop_ratio",
    "ndrugactT",
    "Cl_blin",
    "Vz_blin",
    "KdrugactT_blin",
    "ndrugactT_blin",
    "KmTB_kill_blin",
    "KdrugB_blin",
    "nkill_blin",
    "KmBT_act_blin",
    "kBapop_cll",
    "kBgen_cll",
    "fBkill",
    "fTact",
    "fKmTB_kill",
    "Kptumor",
    "Vtumor",
    "KBptumor",
    "KTrptumor",
    "tumor_on",
    "IL6_tiss_contribution",
    "kBtumorprolif",
    "S",
    "tissue1on",
    "kabs_TDB",
    "fbio_TDB",
    "Bcell_tumor_trafficking_on",
]

function direct_model_layout_matches(state_names::Vector{String}, param_names::Vector{String})
    return state_names == DIRECT_STATE_NAMES && param_names == DIRECT_PARAM_NAMES
end

function apply_initial_rules_direct!(z, u0, pvals)
    @inbounds begin
        z[1:DIRECT_STATE_COUNT] .= u0
        z[(DIRECT_STATE_COUNT + 1):(DIRECT_STATE_COUNT + DIRECT_PARAM_COUNT)] .= pvals
        z[9] = Float64(real((1 - z[195]) * z[153] * z[160]))
        u0[9] = z[9]
        z[7] = Float64(real(z[152] * z[153] * z[161]))
        u0[7] = z[7]
        z[42] = Float64(real(z[165] * z[153] * z[184]))
        u0[42] = z[42]
        z[63] = Float64(real(z[203] * z[153] * z[202]))
        u0[63] = z[63]
        z[102] = Float64(real(z[230] * z[153] * z[228]))
        u0[102] = z[102]
        z[10] = Float64(real((1 - z[196]) * z[156] * z[160]))
        u0[10] = z[10]
        z[2] = Float64(real(z[158] * z[156] * z[161]))
        u0[2] = z[2]
        z[36] = Float64(real(z[183] * z[156] * z[184]))
        u0[36] = z[36]
        z[61] = Float64(real(z[174] * z[207] * z[202]))
        u0[61] = z[61]
        z[62] = Float64(real(z[174] * z[207] * z[208] * z[202]))
        u0[62] = z[62]
        z[104] = Float64(real(z[229] * z[156] * z[228]))
        u0[104] = z[104]
        z[31] = Float64(real(z[179]))
        u0[31] = z[31]
    end
    return nothing
end

function apply_repeated_rules_direct!(z, u, pvals, t)
    @inbounds begin
        z[1:DIRECT_STATE_COUNT] .= u
        z[(DIRECT_STATE_COUNT + 1):(DIRECT_STATE_COUNT + DIRECT_PARAM_COUNT)] .= pvals
        z[11] = Float64(real(PK_v26(z[3], z[148], z[198], z[199], t, z[200], z[209])))
        z[57] = Float64(real(((z[55] / z[187] > 1.0e-5) * z[55]) / z[187]))
        z[85] = Float64(real((z[84] / z[215]) * 1000))
        z[18] = Float64(real(z[9] / z[160]))
        z[15] = Float64(real(z[7] / z[161]))
        z[54] = Float64(real(z[42] / z[184]))
        z[70] = Float64(real(z[63] / z[202]))
        z[111] = Float64(real(z[102] / z[228]))
        z[29] = Float64(real(z[33] / z[160]))
        z[28] = Float64(real(z[34] / z[161]))
        z[48] = Float64(real(z[41] / z[184]))
        z[71] = Float64(real(z[64] / z[202]))
        z[112] = Float64(real(z[105] / z[228]))
        z[19] = Float64(real(z[4] / z[160]))
        z[16] = Float64(real(z[1] / z[161]))
        z[53] = Float64(real(z[46] / z[184]))
        z[72] = Float64(real(z[65] / z[202]))
        z[113] = Float64(real(z[103] / z[228]))
        z[69] = Float64(real(z[62] + z[61]))
        z[20] = Float64(real(z[10] / z[160]))
        z[17] = Float64(real(z[2] / z[161]))
        z[52] = Float64(real(z[36] / z[184]))
        z[81] = Float64(real(z[61] / z[202]))
        z[83] = Float64(real(z[62] / z[202]))
        z[114] = Float64(real(z[104] / z[228]))
        z[21] = Float64(real(z[10] / max(z[9] + z[33], 1)))
        z[22] = Float64(real(z[2] / max(z[7] + z[34], 1)))
        z[51] = Float64(real(z[36] / max(z[42] + z[41], 1)))
        z[82] = Float64(real(z[61] / max(z[63] + z[64], 1)))
        z[116] = Float64(real(z[104] / max(z[102] + z[105], 1)))
        z[23] = Float64(real(z[4] / max(z[10], 1)))
        z[24] = Float64(real(z[1] / max(z[2], 1)))
        z[50] = Float64(real(z[46] / max(z[36], 1)))
        z[80] = Float64(real(z[65] / max(z[61], 1)))
        z[117] = Float64(real(z[103] / max(z[104], 1)))
        z[99] = Float64(real(z[7] + z[34] + z[1]))
        z[100] = Float64(real(z[42] + z[41] + z[46]))
        z[101] = Float64(real(z[63] + z[64] + z[65]))
        z[118] = Float64(real(z[102] + z[105] + z[103]))
        z[32] = Float64(real(((log(2) / z[175]) * z[31] * (z[2] + z[10] + z[36])) / (z[174] * (z[160] + z[158] * z[161] + z[183] * z[184]))))
        z[131] = Float64(real(z[126] + (z[232] * (z[127] * z[161] + z[128] * z[184] + z[129] * z[202] + z[130] * z[228])) / z[160]))
        z[12] = Float64(real(z[159] * z[11]))
        z[43] = Float64(real(z[185] * z[11]))
        z[67] = Float64(real(z[204] * z[205] * z[11]))
        z[110] = Float64(real(z[231] * z[227] * z[11]))
        z[58] = Float64(real(z[57] * z[159]))
        z[59] = Float64(real(z[57] * z[185]))
        z[78] = Float64(real(z[204] * z[57] * z[205]))
        z[121] = Float64(real(z[231] * z[57] * z[227]))
        z[90] = Float64(real(z[159] * z[85]))
        z[91] = Float64(real(z[185] * z[85]))
        z[95] = Float64(real(z[204] * z[205] * z[85]))
        z[120] = Float64(real(z[231] * z[227] * z[85]))
        z[26] = Float64(real(z[18] + z[29] + z[19]))
        z[27] = Float64(real(z[15] + z[28] + z[16]))
        z[49] = Float64(real(z[54] + z[48] + z[53]))
        z[75] = Float64(real(z[70] + z[71] + z[72]))
        z[119] = Float64(real(z[111] + z[112] + z[113]))
        z[73] = Float64(real(z[69] / z[202]))
        z[68] = Float64(real(z[69] / max(z[63] + z[64], 1)))
        z[76] = Float64(real(z[65] / max(z[69], 1)))
        z[97] = Float64(real(z[20] / z[174]))
        z[13] = Float64(real(z[140] * (pow_safe(z[21], z[234]) / (pow_safe(z[141], z[234]) + pow_safe(z[21], z[234]))) * (pow_safe(z[11] * 1000, z[213]) / (pow_safe(z[139], z[213]) + pow_safe(z[11] * 1000, z[213])))))
        z[86] = Float64(real(z[140] * (pow_safe(z[21], z[234]) / (pow_safe(z[221], z[234]) + pow_safe(z[21], z[234]))) * (pow_safe(z[85], z[217]) / (pow_safe(z[85], z[217]) + pow_safe(z[216], z[217])))))
        z[14] = Float64(real(((z[132] * pow_safe(max(0, z[23]), z[164])) / (pow_safe(max(0, z[133]), z[164]) + pow_safe(max(0, z[23]), z[164]))) * ((z[11] * 1000) / (z[134] + z[11] * 1000))))
        z[87] = Float64(real(((z[132] * pow_safe(max(0, z[23]), z[220])) / (pow_safe(max(0, z[218]), z[220]) + pow_safe(max(0, z[23]), z[220]))) * (z[85] / (z[219] + z[85]))))
        z[6] = Float64(real(z[225] * z[140] * (pow_safe(z[22], z[234]) / (pow_safe(z[141], z[234]) + pow_safe(z[22], z[234]))) * (pow_safe(z[12] * 1000, z[213]) / (pow_safe(z[139], z[213]) + pow_safe(z[12] * 1000, z[213])))))
        z[5] = Float64(real(((z[132] * pow_safe(max(0, z[24]), z[164])) / (pow_safe(max(0, z[133] * z[226]), z[164]) + pow_safe(max(0, z[24]), z[164]))) * ((z[12] * 1000) / (z[134] + z[12] * 1000))))
        z[44] = Float64(real(z[225] * z[140] * (pow_safe(z[51], z[234]) / (pow_safe(z[141], z[234]) + pow_safe(z[51], z[234]))) * (pow_safe(z[43] * 1000, z[213]) / (pow_safe(z[139], z[213]) + pow_safe(z[43] * 1000, z[213])))))
        z[45] = Float64(real(((z[132] * pow_safe(max(0, z[50]), z[164])) / (pow_safe(max(0, z[133] * z[226]), z[164]) + pow_safe(max(0, z[50]), z[164]))) * ((z[43] * 1000) / (z[134] + z[43] * 1000))))
        z[66] = Float64(real(z[225] * z[140] * (pow_safe(z[82], z[234]) / (pow_safe(z[141], z[234]) + pow_safe(z[82], z[234]))) * (pow_safe(z[67] * 1000, z[213]) / (pow_safe(z[139], z[213]) + pow_safe(z[67] * 1000, z[213])))))
        z[79] = Float64(real(((z[132] * pow_safe(max(0, z[80]), z[164])) / (pow_safe(max(0, z[133] * z[226]), z[164]) + pow_safe(max(0, z[80]), z[164]))) * ((z[67] * 1000) / (z[134] + z[67] * 1000))))
        z[106] = Float64(real(z[225] * z[140] * (pow_safe(max(z[116], 0), z[234]) / (pow_safe(z[141], z[234]) + pow_safe(max(z[116], 0), z[234]))) * (pow_safe(z[110] * 1000, z[213]) / (pow_safe(z[139], z[213]) + pow_safe(z[110] * 1000, z[213])))))
        z[107] = Float64(real(((z[132] * pow_safe(max(0, z[117]), z[164])) / (pow_safe(max(0, z[133] * z[226]), z[164]) + pow_safe(max(0, z[117]), z[164]))) * ((z[110] * 1000) / (z[134] + z[110] * 1000))))
        z[88] = Float64(real(z[225] * z[140] * (pow_safe(z[22], z[234]) / (pow_safe(z[221], z[234]) + pow_safe(z[22], z[234]))) * (pow_safe(z[90], z[217]) / (pow_safe(z[90], z[217]) + pow_safe(z[216], z[217])))))
        z[89] = Float64(real(((z[132] * pow_safe(max(0, z[24]), z[220])) / (pow_safe(max(0, z[218] * z[226]), z[220]) + pow_safe(max(0, z[24]), z[220]))) * (z[90] / (z[219] + z[90]))))
        z[93] = Float64(real(z[225] * z[140] * (pow_safe(z[51], z[234]) / (pow_safe(z[221], z[234]) + pow_safe(z[51], z[234]))) * (pow_safe(z[91], z[217]) / (pow_safe(z[91], z[217]) + pow_safe(z[216], z[217])))))
        z[92] = Float64(real(((z[132] * pow_safe(max(0, z[50]), z[220])) / (pow_safe(max(0, z[218] * z[226]), z[220]) + pow_safe(max(0, z[50]), z[220]))) * (z[91] / (z[219] + z[91]))))
        z[108] = Float64(real(z[225] * z[140] * (pow_safe(max(z[116], 0), z[234]) / (pow_safe(z[221], z[234]) + pow_safe(max(z[116], 0), z[234]))) * (pow_safe(z[120], z[217]) / (pow_safe(z[120], z[217]) + pow_safe(z[216], z[217])))))
        z[109] = Float64(real(((z[132] * pow_safe(max(0, z[117]), z[220])) / (pow_safe(max(0, z[218] * z[226]), z[220]) + pow_safe(max(0, z[117]), z[220]))) * (z[120] / (z[219] + z[120]))))
        z[25] = Float64(real(z[19] / max(z[26], 1)))
        z[30] = Float64(real(z[16] / max(z[27], 1)))
        z[122] = Float64(real(z[17] / z[27]))
        z[47] = Float64(real(z[53] / max(z[49], 1)))
        z[123] = Float64(real(z[52] / z[49]))
        z[74] = Float64(real(z[72] / max(z[75], 1)))
        z[115] = Float64(real(z[113] / max(z[119], 1)))
        z[96] = Float64(real(z[225] * z[140] * (pow_safe(z[68], z[234]) / (pow_safe(z[221], z[234]) + pow_safe(z[68], z[234]))) * (pow_safe(z[95], z[217]) / (pow_safe(z[95], z[217]) + pow_safe(z[216], z[217])))))
        z[94] = Float64(real(((z[132] * pow_safe(max(0, z[76]), z[220])) / (pow_safe(max(0, z[218] * z[226]), z[220]) + pow_safe(max(0, z[76]), z[220]))) * (z[95] / (z[219] + z[95]))))
    end
    return z
end

function rhs_direct!(du, u, ctx, t)
    z = ctx.z
    pvals = ctx.pvals
    apply_repeated_rules_direct!(z, u, pvals, t)
    @inbounds begin
        fill!(du, 0.0)
        for inf in ctx.infusions
            if t >= inf.t_start && t <= inf.t_end
                du[inf.target_idx] += inf.rate
            end
        end
        rate_1 = Float64(real(z[136] * z[158] * z[174] * z[161] * pow_safe(max(0, 1 - z[2] / (z[158] * z[174] * z[161])), 1)))
        du[2] += 1.0 * rate_1
        rate_2 = Float64(real(z[143] * ((z[6] + z[88]) * z[7] - z[144] * z[1])))
        du[1] += 1.0 * rate_2
        du[7] += -1.0 * rate_2
        rate_3 = Float64(real(0 * z[155] * z[136] * z[2] * z[222]))
        du[2] += -1.0 * rate_3
        rate_4 = Float64(real(((z[146] + (z[166] / (z[167] + z[3] / z[148])) / z[168]) / z[148]) * z[3]))
        du[3] += -1.0 * rate_4
        rate_5 = Float64(real(z[147] * (z[3] / z[148] - z[8] / z[149])))
        du[3] += -1.0 * rate_5
        du[8] += 1.0 * rate_5
        rate_6 = Float64(real(0))
        du[7] += 1.0 * rate_6
        rate_7 = Float64(real(0))
        rate_8 = Float64(real(z[235] * z[197] * z[160] * ((((1 + z[182] * z[35] + z[201] * z[60]) * z[9]) / z[160]) * z[152] - z[7] / z[161])))
        du[7] += 1.0 * rate_8
        du[9] += -1.0 * rate_8
        rate_9 = Float64(real(z[150] * z[160] * z[163] * (z[194] + pow_safe(max(0, 1 - z[26] / z[163]), 2))))
        du[9] += 1.0 * rate_9
        rate_10 = Float64(real(z[235] * z[157] * z[142] * z[160] * max(0, (z[10] / z[160]) * z[158] - z[2] / z[161])))
        du[2] += 1.0 * rate_10
        du[10] += -1.0 * rate_10
        rate_11 = Float64(real(z[143] * ((z[13] + z[86]) * z[9] - z[144] * z[4])))
        du[4] += 1.0 * rate_11
        du[9] += -1.0 * rate_11
        rate_12 = Float64(real(z[204] * z[211] * z[135] * z[10] * z[222]))
        du[10] += 1.0 * rate_12
        rate_13 = Float64(real(z[135] * z[10] * z[222]))
        du[10] += -1.0 * rate_13
        rate_14 = Float64(real(0))
        rate_15 = Float64(real(0))
        rate_16 = Float64(real(z[138] * (z[14] + z[87] + z[57] / (z[193] / 1000 + z[57])) * z[10]))
        du[10] += -1.0 * rate_16
        rate_17 = Float64(real(z[224] * z[138] * (z[5] + z[89] + z[58] / (z[193] / 1000 + z[58])) * z[2]))
        du[2] += -1.0 * rate_17
        rate_18 = Float64(real(z[150] * z[194] * z[9] + z[162] * z[151] * max(0, z[9] - z[163] * z[160])))
        du[9] += -1.0 * rate_18
        rate_19 = Float64(real(z[235] * z[142] * z[160] * ((((1 + z[171] * (z[182] * z[35] + z[201] * z[60])) * z[4]) / z[160]) * z[152] * z[145] - z[1] / z[161])))
        du[1] += 1.0 * rate_19
        du[4] += -1.0 * rate_19
        rate_20 = Float64(real(z[151] * (z[4] + (z[172] * pow_safe(z[4], 2)) / (z[160] * z[163]))))
        du[4] += -1.0 * rate_20
        rate_21 = Float64(real((log(2) / ((z[175] / 24) / 60)) * ((z[180] * z[31]) / z[179] + ((1 - z[180]) * (z[2] + z[10] + z[36])) / (z[174] * (z[160] + z[158] * z[161] + z[183] * z[184])))))
        du[31] += -1.0 * rate_21
        rate_22 = Float64(real((log(2) / ((z[175] / 24) / 60)) * (z[180] + ((1 - z[180]) * z[156]) / z[174])))
        du[31] += 1.0 * rate_22
        rate_23 = Float64(real(z[169] * z[177] * z[33]))
        du[9] += 1.0 * rate_23
        du[33] += -1.0 * rate_23
        rate_24 = Float64(real(z[169] * z[143] * ((z[13] + z[86]) * z[33] - z[144] * z[4])))
        du[4] += 1.0 * rate_24
        du[33] += -1.0 * rate_24
        rate_25 = Float64(real(z[178] * z[151] * z[33]))
        du[33] += -1.0 * rate_25
        rate_26 = Float64(real(z[178] * z[151] * z[34]))
        du[34] += -1.0 * rate_26
        rate_27 = Float64(real(z[169] * z[143] * ((z[6] + z[88]) * z[34] - z[144] * z[1])))
        du[1] += 1.0 * rate_27
        du[34] += -1.0 * rate_27
        rate_28 = Float64(real(z[169] * z[177] * z[34]))
        du[7] += 1.0 * rate_28
        du[34] += -1.0 * rate_28
        rate_29 = Float64(real(z[235] * z[169] * z[197] * z[160] * ((((1 + z[182] * z[35] + z[201] * z[60]) * z[33]) / z[160]) * z[152] - z[34] / z[161])))
        du[33] += -1.0 * rate_29
        du[34] += 1.0 * rate_29
        rate_30 = Float64(real(z[169] * z[154] * z[137] * z[4]))
        du[33] += 1.0 * rate_30
        rate_31 = Float64(real(z[169] * z[154] * z[137] * z[1]))
        du[34] += 1.0 * rate_31
        rate_32 = Float64(real((log(2) / z[181]) * z[35]))
        du[35] += -1.0 * rate_32
        rate_33 = Float64(real(z[169] * z[154] * z[137] * z[46]))
        du[41] += 1.0 * rate_33
        rate_34 = Float64(real(z[186] * z[169] * z[197] * z[160] * ((((1 + z[182] * z[35] + z[201] * z[60]) * z[33]) / z[160]) * z[165] - z[41] / z[184])))
        du[33] += -1.0 * rate_34
        du[41] += 1.0 * rate_34
        rate_35 = Float64(real(z[169] * z[177] * z[41]))
        du[41] += -1.0 * rate_35
        du[42] += 1.0 * rate_35
        rate_36 = Float64(real(z[169] * z[143] * ((z[44] + z[93]) * z[41] - z[144] * z[46])))
        du[41] += -1.0 * rate_36
        du[46] += 1.0 * rate_36
        rate_37 = Float64(real(z[178] * z[151] * z[41]))
        du[41] += -1.0 * rate_37
        rate_38 = Float64(real(z[186] * z[142] * z[160] * ((((1 + z[171] * (z[182] * z[35] + z[201] * z[60])) * z[4]) / z[160]) * z[165] * z[145] - z[46] / z[184])))
        du[4] += -1.0 * rate_38
        du[46] += 1.0 * rate_38
        rate_39 = Float64(real(z[224] * z[138] * (z[45] + z[92] + z[59] / (z[193] / 1000 + z[59])) * z[36]))
        du[36] += -1.0 * rate_39
        rate_40 = Float64(real(z[186] * z[157] * z[142] * z[160] * max(0, (z[10] / z[160]) * z[183] - z[36] / z[184])))
        du[10] += -1.0 * rate_40
        du[36] += 1.0 * rate_40
        rate_41 = Float64(real(z[186] * z[197] * z[160] * ((((1 + z[182] * z[35] + z[201] * z[60]) * z[9]) / z[160]) * z[165] - z[42] / z[184])))
        du[9] += -1.0 * rate_41
        du[42] += 1.0 * rate_41
        rate_42 = Float64(real(0))
        rate_43 = Float64(real(0))
        rate_44 = Float64(real(0 * z[155] * z[136] * z[36] * z[222]))
        du[36] += -1.0 * rate_44
        rate_45 = Float64(real(z[143] * ((z[44] + z[93]) * z[42] - z[144] * z[46])))
        du[42] += -1.0 * rate_45
        du[46] += 1.0 * rate_45
        rate_46 = Float64(real(z[136] * z[183] * z[174] * z[184] * pow_safe(max(0, 1 - z[36] / (z[183] * z[174] * z[184])), 1)))
        du[36] += 1.0 * rate_46
        rate_47 = Float64(real(z[190] * (z[55] / z[187] - z[56] / z[188])))
        du[55] += -1.0 * rate_47
        du[56] += 1.0 * rate_47
        rate_48 = Float64(real(((z[189] + (z[191] / (z[192] + z[55] / z[187])) / z[168]) / z[187]) * z[55]))
        du[55] += -1.0 * rate_48
        rate_49 = Float64(real((log(2) / z[181]) * z[60]))
        du[60] += -1.0 * rate_49
        rate_50 = Float64(real(z[174] * z[207] * z[208] * z[202] * (z[135] + z[135] / z[212]) * (1 + 5 * pow_safe(max(0, 1 - (z[10] + z[2] + z[36]) / (z[174] * z[160] + z[158] * z[156] * z[161] + z[183] * z[156] * z[184])), 1))))
        du[62] += 1.0 * rate_50
        rate_51 = Float64(real((z[135] / z[212]) * z[62]))
        du[61] += 1.0 * rate_51
        du[62] += -1.0 * rate_51
        rate_52 = Float64(real(z[204] * z[211] * z[206] * z[160] * max(0, z[61] / z[202] - (z[10] / z[160]) * z[207])))
        du[10] += 1.0 * rate_52
        du[61] += -1.0 * rate_52
        rate_53 = Float64(real(z[135] * z[62] * z[222]))
        du[62] += -1.0 * rate_53
        rate_54 = Float64(real(z[135] * z[61] * z[222]))
        du[61] += -1.0 * rate_54
        rate_55 = Float64(real(z[204] * z[197] * z[160] * ((((1 + z[182] * z[35] + z[201] * z[60]) * z[9]) / z[160]) * z[203] - z[63] / z[202])))
        du[9] += -1.0 * rate_55
        du[63] += 1.0 * rate_55
        rate_56 = Float64(real(z[224] * z[138] * (z[79] + z[94] + z[78] / (z[193] / 1000 + z[78])) * z[61]))
        du[61] += -1.0 * rate_56
        rate_57 = Float64(real(z[204] * z[169] * z[197] * z[160] * ((((1 + z[182] * z[35] + z[201] * z[60]) * z[33]) / z[160]) * z[203] - z[64] / z[202])))
        du[33] += -1.0 * rate_57
        du[64] += 1.0 * rate_57
        rate_58 = Float64(real(z[204] * z[142] * z[160] * ((((1 + z[171] * (z[182] * z[35] + z[201] * z[60])) * z[4]) / z[160]) * z[203] * z[145] - z[65] / z[202])))
        du[4] += -1.0 * rate_58
        du[65] += 1.0 * rate_58
        rate_59 = Float64(real(z[178] * z[151] * z[64]))
        du[64] += -1.0 * rate_59
        rate_60 = Float64(real(z[169] * z[143] * ((z[66] + z[96]) * z[64] - z[144] * z[65])))
        du[64] += -1.0 * rate_60
        du[65] += 1.0 * rate_60
        rate_61 = Float64(real(z[169] * z[177] * z[64]))
        du[63] += 1.0 * rate_61
        du[64] += -1.0 * rate_61
        rate_62 = Float64(real(z[143] * ((z[66] + z[96]) * z[63] - z[144] * z[65])))
        du[63] += -1.0 * rate_62
        du[65] += 1.0 * rate_62
        rate_63 = Float64(real(z[169] * z[154] * z[137] * z[65]))
        du[64] += 1.0 * rate_63
        rate_64 = Float64(real(0))
        rate_65 = Float64(real(z[224] * z[138] * z[94] * z[62]))
        du[62] += -1.0 * rate_65
        rate_66 = Float64(real(0))
        rate_67 = Float64(real(0))
        rate_68 = Float64(real(z[210] * z[162] * z[151] * (z[163] * z[161] * z[152]) * max(0, z[7] / (z[163] * z[161] * z[152]) - 1)))
        du[7] += -1.0 * rate_68
        rate_69 = Float64(real(z[210] * z[162] * z[151] * (z[163] * z[184] * z[165]) * max(0, z[42] / (z[163] * z[184] * z[165]) - 1)))
        du[42] += -1.0 * rate_69
        rate_70 = Float64(real(z[210] * z[151] * (z[1] + (z[172] * pow_safe(z[1], 2)) / (z[152] * z[161] * z[163]))))
        du[1] += -1.0 * rate_70
        rate_71 = Float64(real(z[210] * z[151] * (z[46] + (z[172] * pow_safe(z[46], 2)) / (z[184] * z[165] * z[163]))))
        du[46] += -1.0 * rate_71
        rate_72 = Float64(real(z[210] * z[162] * z[151] * (z[163] * z[202] * z[203]) * max(0, z[63] / (z[163] * z[202] * z[203]) - 1)))
        du[63] += -1.0 * rate_72
        rate_73 = Float64(real(z[210] * z[151] * (z[65] + (z[172] * pow_safe(z[65], 2)) / (z[202] * z[203] * z[163]))))
        du[65] += -1.0 * rate_73
        rate_74 = Float64(real(z[204] * z[211] * z[136] * z[222] * z[174] * z[160] * pow_safe(max(0, 1 - z[10] / (z[174] * z[160])), 1)))
        du[10] += 1.0 * rate_74
        rate_75 = Float64(real((z[214] * z[84]) / z[215]))
        du[84] += -1.0 * rate_75
        rate_76 = Float64(real(0))
        rate_77 = Float64(real(0))
        rate_78 = Float64(real(0))
        rate_79 = Float64(real(0))
        rate_80 = Float64(real(0))
        rate_81 = Float64(real(0))
        rate_82 = Float64(real(0))
        rate_83 = Float64(real(0))
        rate_84 = Float64(real(z[136] * z[207] * z[174] * z[202] * pow_safe(max(0, 1 - z[61] / (z[207] * z[174] * z[202])), 1)))
        du[61] += 1.0 * rate_84
        rate_85 = Float64(real(z[136] * z[207] * z[174] * z[208] * z[202] * pow_safe(max(0, 1 - z[62] / (z[207] * z[174] * z[208] * z[202])), 1)))
        du[62] += 1.0 * rate_85
        rate_86 = Float64(real(z[169] * z[177] * z[105]))
        du[102] += 1.0 * rate_86
        du[105] += -1.0 * rate_86
        rate_87 = Float64(real(z[169] * z[143] * ((z[106] + z[108]) * z[105] - z[144] * z[103])))
        du[103] += 1.0 * rate_87
        du[105] += -1.0 * rate_87
        rate_88 = Float64(real(z[169] * z[154] * z[137] * z[103]))
        du[105] += 1.0 * rate_88
        rate_89 = Float64(real(z[143] * ((z[106] + z[108]) * z[102] - z[144] * z[103])))
        du[102] += -1.0 * rate_89
        du[103] += 1.0 * rate_89
        rate_90 = Float64(real(0))
        rate_91 = Float64(real(0))
        rate_92 = Float64(real(0))
        rate_93 = Float64(real(0))
        rate_94 = Float64(real(z[233] * z[104] + 0 * z[136] * z[229] * z[174] * z[228] * pow_safe(max(0, 1 - z[104] / (z[229] * z[174] * z[228])), 1)))
        du[104] += 1.0 * rate_94
        rate_95 = Float64(real(0 * z[155] * z[136] * z[104]))
        du[104] += -1.0 * rate_95
        rate_96 = Float64(real(z[224] * z[138] * (z[107] + z[109] + z[121] / (z[193] / 1000 + z[121])) * z[104]))
        du[104] += -1.0 * rate_96
        rate_97 = Float64(real(z[210] * z[151] * (z[103] + (z[172] * pow_safe(z[103], 2)) / (z[228] * z[230] * z[163]))))
        du[103] += -1.0 * rate_97
        rate_98 = Float64(real(z[238] * z[231] * z[157] * z[142] * z[160] * max(0, (z[10] / z[160]) * z[229] - z[104] / z[228])))
        du[10] += -1.0 * rate_98
        du[104] += 1.0 * rate_98
        rate_99 = Float64(real(z[231] * z[197] * z[160] * ((((1 + z[182] * z[35] + z[201] * z[60]) * z[9]) / z[160]) * z[230] - z[102] / z[228])))
        du[9] += -1.0 * rate_99
        du[102] += 1.0 * rate_99
        rate_100 = Float64(real(z[231] * z[142] * z[160] * ((((1 + z[171] * (z[182] * z[35] + z[201] * z[60])) * z[4]) / z[160]) * z[230] * z[145] - z[103] / z[228])))
        du[4] += -1.0 * rate_100
        du[103] += 1.0 * rate_100
        rate_101 = Float64(real(z[231] * z[169] * z[197] * z[160] * ((((1 + z[182] * z[35] + z[201] * z[60]) * z[33]) / z[160]) * z[230] - z[105] / z[228])))
        du[33] += -1.0 * rate_101
        du[105] += 1.0 * rate_101
        rate_102 = Float64(real(z[178] * z[151] * z[105]))
        du[105] += -1.0 * rate_102
        rate_103 = Float64(real(z[210] * z[162] * z[151] * (z[163] * z[228] * z[230]) * max(0, z[102] / (z[163] * z[228] * z[230]) - 1)))
        du[102] += -1.0 * rate_103
        rate_104 = Float64(real(z[236] * z[237] * z[124]))
        du[3] += 1.0 * rate_104
        du[124] += -1.0 * rate_104
        rate_105 = Float64(real(z[236] * (1 - z[237]) * z[124]))
        du[124] += -1.0 * rate_105
        rate_106 = Float64(real(z[11]))
        du[125] += 1.0 * rate_106
        rate_107 = Float64(real(((z[170] * z[4]) / z[160]) * (z[20] / z[174]) * (pow_safe(z[11] * 1000, z[213]) / (pow_safe(z[139], z[213]) + pow_safe(z[11] * 1000, z[213])) + pow_safe(z[85], z[217]) / (pow_safe(z[85], z[217]) + pow_safe(z[216], z[217])))))
        du[126] += 1.0 * rate_107
        rate_108 = Float64(real((log(2) / ((z[176] / 60) / 24)) * z[126]))
        du[126] += -1.0 * rate_108
        rate_109 = Float64(real(((z[170] * z[1]) / z[161]) * (z[17] / (z[174] * z[158])) * (pow_safe(z[12] * 1000, z[213]) / (pow_safe(z[139], z[213]) + pow_safe(z[12] * 1000, z[213])) + pow_safe(z[90], z[217]) / (pow_safe(z[90], z[217]) + pow_safe(z[216], z[217])))))
        du[127] += 1.0 * rate_109
        rate_110 = Float64(real((log(2) / ((z[176] / 60) / 24)) * z[127]))
        du[127] += -1.0 * rate_110
        rate_111 = Float64(real(((z[170] * z[46]) / z[184]) * (z[52] / (z[174] * z[183])) * (pow_safe(z[43] * 1000, z[213]) / (pow_safe(z[139], z[213]) + pow_safe(z[43] * 1000, z[213])) + pow_safe(z[91], z[217]) / (pow_safe(z[91], z[217]) + pow_safe(z[216], z[217])))))
        du[128] += 1.0 * rate_111
        rate_112 = Float64(real((log(2) / ((z[176] / 60) / 24)) * z[128]))
        du[128] += -1.0 * rate_112
        rate_113 = Float64(real(((z[170] * z[65]) / z[202]) * ((z[81] / (z[174] * z[207])) * (pow_safe(z[67] * 1000, z[213]) / (pow_safe(z[139], z[213]) + pow_safe(z[67] * 1000, z[213]))) + (z[73] / (z[174] * z[207] * (1 + z[208]))) * (pow_safe(z[95], z[217]) / (pow_safe(z[95], z[217]) + pow_safe(z[216], z[217]))))))
        du[129] += 1.0 * rate_113
        rate_114 = Float64(real((log(2) / ((z[176] / 60) / 24)) * z[129]))
        du[129] += -1.0 * rate_114
        rate_115 = Float64(real(((z[170] * z[103]) / z[228]) * (z[114] / (z[174] * z[229])) * (pow_safe(z[110] * 1000, z[213]) / (pow_safe(z[139], z[213]) + pow_safe(z[110] * 1000, z[213])) + pow_safe(z[120], z[217]) / (pow_safe(z[120], z[217]) + pow_safe(z[216], z[217])))))
        du[130] += 1.0 * rate_115
        rate_116 = Float64(real((log(2) / ((z[176] / 60) / 24)) * z[130]))
        du[130] += -1.0 * rate_116
    end
    return nothing
end
