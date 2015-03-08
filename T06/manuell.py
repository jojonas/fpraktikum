from collections import OrderedDict

kalib = OrderedDict()

kalib["G20_Kalib_Cu.mca"] = {
	"title": "Kalibrationstarget Kupfer (Cu, Z = 29)",
	"element": "Cu",
	"peaks": [
		(600, 10, r'K_\alpha', 8.048, 100),
		(665, 10, r'K_\beta', 8.903, 6.21),
	]
}

kalib["G20_Kalib_Rb.mca"] = {
	"title": "Kalibrationstarget Rubidium (Rb, Z = 37)",
	"element": "Rb",
	"peaks": [
		(1000, 10, r'K_\alpha', 13.395, 100),
		(1115, 10, r'K_\beta', 14.961, 14.52),
	]
}

kalib["G20_Kalib_Mo.mca"] = {
	"title": "Kalibrationstarget Molybdän (Mo, Z = 42)",
	"element": "Mo",
	"peaks": [
		(1300, 15, r'K_\alpha', 17.479, 100),
		(1465, 10, r'K_\beta', 19.608, 15.86),
	]
}

kalib["G20_Kalib_Ag.mca"] = {
	"title": "Kalibrationstarget Silber (Ag, Z = 47)",
	"element": "Ag",
	"peaks": [
		(1650, 10, r'K_\alpha', 22.163, 100),
		(1865, 10, r'K_\beta', 24.942, 16.95),
	]
}

kalib["G20_Kalib_Ba.mca"] = {
 	"title": "Kalibrationstarget Barium (Ba, Z = 56)",
	"element": "Ba",
	"peaks": [
		# (2410, 10, 32.06), # K-alpha (Liste)
		# (2720, 20, 36.55), # K-beta (Liste)
		#(2380, 10, 31.82), # KL2
		(2410, 10, r'K_{\alpha 1}', 32.194, 100),
		(2720, 20, r'K_{\beta 1}', 36.378, 18.53),
		(2790, 10, r'K_{\beta 2}', 37.257, 3.87),
	]
}

kalib["G20_Kalib_Tb.mca"] = {
	"title": "Kalibrationstarget Terbium (Tb, Z = 65)",
	"element": "Tb",
	"peaks": [
		# L alpha 2: 11.32 % 6238eV
		#(470, 10, r'L_{\alpha 1}', 6.272), # 100%
#		(520, 10, r'L_{\alpha 5}', 7.01),
		#(550, 10, r'L_{\beta 4/5}^*', 7.37), # 70.90
		#(550, 10, r'L_{\beta 4}^*', 6.940), # 70.90 ???????
		#(605, 5, r'L_{\gamma 1}', 8.102), #17.46%
		# L-gamma 2 8398eV 16.79%
		# L gamma3 8423eV 24-68%
		(3275, 20, r'K_{\alpha 1}', 43.744, 100), #55.74),
		(3335, 15, r'K_{\alpha 2}', 44.482, 179.4), #100),
		#(3770, 20, r'K_{\beta 3}^*', 50.229), #10.22%
		(3770, 20, r'K_{\beta 1}', 50.382, 35.47), #19.77),
		(3875, 15, r'K_{\beta 2}', 51.720, 7.8) #4.35),
		# Lbeta 2 7366.7 # 16.87%
	]
}

data = OrderedDict()

data["G20_Unbekannt1_2cent.mca"] =	{
	"title": "2-Cent-Münze",
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
