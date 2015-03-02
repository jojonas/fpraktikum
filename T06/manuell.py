from collections import OrderedDict

kalib = OrderedDict()

kalib["G20_Kalib_Cu.mca"] = {
	"title": "Kalibrationstarget Kupfer (Cu, Z = 29)",
	"peaks": [
		(600, 10, 8.04),
		(665, 10, 8.91)
	]
}

kalib["G20_Kalib_Rb.mca"] = {
	"title": "Kalibrationstarget Rubidium (Rb, T = 37)",
	"peaks": [
		(1000, 10, 13.37),
		(1115, 10, 14.97)
	]
}

kalib["G20_Kalib_Mo.mca"] = {
	"title": "Kalibrationstarget Molybdän (Mo, Z = 42)",
	"peaks": [
		(1300, 15, 17.44),
		(1465, 10, 19.63)
	]
}

kalib["G20_Kalib_Ag.mca"] = {
	"title": "Kalibrationstarget Silber (Ag, Z = 47)",
	"peaks": [
		(1650, 10, 22.1),
		(1865, 10, 24.99)
	]
}

kalib["G20_Kalib_Ba.mca"] = {
 	"title": "Kalibrationstarget Barium (Ba, Z = 56)",
	"peaks": [
		# (2410, 10, 32.06), # K-alpha (Liste)
		# (2720, 20, 36.55), # K-beta (Liste)
		#(2380, 10, 31.82), # KL2
		(2410, 10, 32.19), # KL3
		(2720, 20, 36.38), # KM3
		(2790, 10, 37.26), # KN2/3
	]
}

kalib["G20_Kalib_Tb.mca"] = {
	"title": "Kalibrationstarget Terbium (Tb, Z = 65)",
	"peaks": [
		# (3275, 20, 44.23), # K-alpha (Liste)
		# (3770, 20, 50.65), # K-beta (Liste)
		(470, 10, 6.28), # L2M1
		(520, 10, 7.01), # L2M5
		(550, 10, 7.37), # L3N4/5
		(605, 5, 8.10),  # L2N4
		(3275, 20, 43.74), # KL2
		(3335, 15, 44.48), # KL3
		(3770, 20, 50.22), # KM3 or KM2
		(3875, 15, 51.72) # KN2/3
	]
}

data = OrderedDict()

data["G20_Unbekannt1_2cent.mca"] =	{
	"title": "Unbekannte Probe: 2-Cent-Stück",
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
	"title": "Unbekannte Probe: Kronkorken",
	"peaks": [
		(480, 10),
		(530, 10),
	]
}

#data["data/G20_Unbekannt4_bleistiftminen.mca"] = {},

data["G20_Unbekannt5_metallplatte.mca"] = {
	"title": "Unbekannte Probe: Metallplatte",
	"peaks" : [
		(600, 10),
		(640, 10)
	]
}
