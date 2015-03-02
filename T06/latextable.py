import codecs

class LatexTable:
	def __init__(self, filename, header=None):
		self.file = codecs.open(filename, 'w', 'utf-8')
		if header:
			self.header(header)

	def close(self):
		if self.file:
			self.file.close()
			self.file = None

	def header(self, *columns, count=None, align=None, borders=True):
		if count is None:
			count = len(columns)

		if align is None:
			align = 'c'
		if not isinstance(align, (list, tuple)):
			align = tuple([align] * count)

		self.file.write(r'\begin{tabular}')
		separator = " | " if borders else " "
		self.file.write(r'{' + separator + separator.join(align) + separator + '}\n')
		self.hline()
		if columns:
			self.row(*columns)
			self.hline(2)

	def footer(self):
		self.hline()
		self.file.write(r'\end{tabular}')
		self.file.write("\n")

	def hline(self, count=1):
		for _ in range(count):
			self.file.write(r'\hline')
		self.file.write("\n")

	def row(self, *values):
		values = (str(v) for v in values)
		separator = " & "
		self.file.write(separator.join(values))
		self.file.write(" \\\\ \n")

	def __enter__(self):
		return self

	def __exit__(self, type, value, traceback):
		self.footer()
		self.close()

	def __del__(self):
		self.close()

if __name__=="__main__":
	import os
	PATH = os.path.dirname(os.path.realpath(__file__))

	with LatexTable(PATH + "/test.tex") as table:
		table.header("Foo", "Bar")
		table.row(1,2)
