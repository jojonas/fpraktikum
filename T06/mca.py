import numpy as np

class McaFile:
	def __init__(self, filename):
		mode = None
		self.filename = filename
		self.data = []
		self.meta = {}
		with open(filename, 'r') as file:
			for line in file:
				line = line.strip()
				if line == "<<DATA>>":
					mode = 'data'
				elif line == "<<PMCA SPECTRUM>>":
					mode = 'pmca'
				elif line in ("<<DPP CONFIGURATION>>", "<<DPP STATUS>>"):
					mode = 'dpp'
				elif line == "<<END>>":
					mode = None
				else:
					if mode == 'data':
						self.data.append(int(line))
					elif mode == 'pmca':
						name, _, value = line.rpartition(" - ")
						self.meta[name] = value
					elif mode == 'dpp':
						name, _, value = line.rpartition(": ")
						self.meta[name] = value

		self.data = np.array(self.data)

	def __len__(self):
		return len(self.data)
