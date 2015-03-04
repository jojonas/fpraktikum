from collections import OrderedDict

kalib = OrderedDict()

kalib["G20_Kalib_Cu.mca"] = {
	"title": "Kalibrationstarget Kupfer (Cu, Z = 29)",
	"element": "Cu",
	"peaks": [
		(600, 10, r'K_\alpha', 8.04),
		(665, 10, r'K_\beta', 8.91)
	]
}

kalib["G20_Kalib_Rb.mca"] = {
	"title": "Kalibrationstarget Rubidium (Rb, Z = 37)",
	"element": "Rb",
	"peaks": [
		(1000, 10, r'K_\alpha', 13.37),
		(1115, 10, r'K_\beta', 14.97)
	]
}

kalib["G20_Kalib_Mo.mca"] = {
	"title": "Kalibrationstarget Molybdän (Mo, Z = 42)",
	"element": "Mo",
	"peaks": [
		(1300, 15, r'K_\alpha', 17.44),
		(1465, 10, r'K_\beta', 19.63)
	]
}

kalib["G20_Kalib_Ag.mca"] = {
	"title": "Kalibrationstarget Silber (Ag, Z = 47)",
	"element": "Ag",
	"peaks": [
		(1650, 10, r'K_\alpha', 22.1),
		(1865, 10, r'K_\beta', 24.99)
	]
}

kalib["G20_Kalib_Ba.mca"] = {
 	"title": "Kalibrationstarget Barium (Ba, Z = 56)",
	"element": "Ba",
	"peaks": [
		# (2410, 10, 32.06), # K-alpha (Liste)
		# (2720, 20, 36.55), # K-beta (Liste)
		#(2380, 10, 31.82), # KL2
		(2410, 10, r'K_{\alpha 3}', 32.19), # KL3 = K-alpha_3
		(2720, 20, r'K_{\beta 3}', 36.38), # KM3 = K-beta_3
		(2790, 10, r'K_{\gamma 2/3}^*', 37.26), # KN2/3 = K-gamma_2 oder K-gamma_3
	]
}

kalib["G20_Kalib_Tb.mca"] = {
	"title": "Kalibrationstarget Terbium (Tb, Z = 65)",
	"element": "Tb",
	"peaks": [
		# (3275, 20, 44.23), # K-alpha (Liste)
		# (3770, 20, 50.65), # K-beta (Liste)
		(470, 10, r'L_{\alpha 1}', 6.28), # L2M1 = L-alpha_1
		(520, 10, r'L_{\alpha 5}', 7.01), # L2M5 = L-alpha_5
		(550, 10, r'L_{\beta 4/5}^*', 7.37), # L3N4/5 = L-beta_4 oder L-beta_5
		(605, 5, r'L_{\gamma 4}', 8.10),  # L2N4 = L-gamma_4
		(3275, 20, r'K_{\alpha 2}', 43.74), # KL2 = K-alpha_2
		(3335, 15, r'K_{\alpha 3}', 44.48), # KL3 = K-alpha_3
		(3770, 20, r'K_{\beta 2/3}^*', 50.22), # KM2/3 = K-beta_2 oder K-beta_3
		(3875, 15, r'K_{\gamma 2/3}^*', 51.72) # KN2/3 = K-gamma_2 oder K-gamma_3
	]
}

data = OrderedDict()

data["G20_Unbekannt1_2cent.mca"] =	{
	"title": "2-Cent-Stück",
	"peaks": [
		(600, 10),
		(665, 10),
	]
}

# data["G20_Unbekannt2_10er.mca"] = {
# 		"peaks": [
# 			(3750, 100)
# 		]
# }

data["G20_Unbekannt3_kronkorken.mca"] = {
	"title": "Kronkorken",
	"peaks": [
		(480, 10),
		(530, 10),
	]
}

#data["data/G20_Unbekannt4_bleistiftminen.mca"] = {},

data["G20_Unbekannt5_metallplatte.mca"] = {
	"title": "Metallplatte",
	"peaks" : [
		(600, 10),
		(640, 10)
	]
}
