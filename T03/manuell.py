from collections import OrderedDict

energy_spectra = OrderedDict()

# energy sources: https://www-nds.iaea.org/xgamma_standards/genergies1.htm

energy_spectra["EnergiespektrumNa.tka"] = {
	"title": "Natrium (Na)",
	"element": "Na",
	"activity": (106.5, 0.8),
	"peaks": [
		(5400, 220, 510.998910, 0.000001),
		(13060, 400, 1274.537, 0.003),
	]
}

energy_spectra["EnergiespektrumCo.tka"] = {
	"title": "Cobalt (Co)",
	"element": "Co",
	"activity": (662, 3),
	"peaks": [
		#(890, 20),
		#(2420, 200),
		(12085, 280, 1173.228, 0.003),
		(13675, 350, 1332.492, 0.004),
	]
}

energy_spectra["EnergiespektrumCs.tka"] = {
	"title": "Cäsium (Cs)",
	"element": "Cs",
	"activity": (36.2, 0.3),
	"peaks": [
		#(388, 40),
		(6930, 220, 661.657, 0.003),
	]
}

energy_spectra["EnergiespektrumEu.tka"] = {
	"title": "Europium (Eu)",
	"element": "Eu",
	"activity": (29.6, 0.2),
	"peaks": [
		#(130, 20),
		#(485, 50),
		(1390, 70, 121.7817, 0.0003),
		(2650, 100, 244.6974, 0.0008),
		(3700, 150, 344.2785, 0.0012),
		(8100, 250, 778.9045, 0.0024),
		(9950, 200, 964.072, 0.018),
		#(11300, 400, 1089.737, 0.005),
		(14360, 320, 1408.013, 0.003),
	]
}
