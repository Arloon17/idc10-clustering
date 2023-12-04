import matplotlib.pyplot as plt

black = [0,0,0]
white = [255,255,255]
gray = [100,100,100]

WINDOW_WIDTH = 1000

def produce_line(length, color):
	line = []
	for _ in range(length):
		line.append(color)
	return line

def color_line(line, color, position):
	line[position] = color
	return line

def color_strip(strip, color, start, end):
	new_strip = []
	for i in range(len(strip)):
		if i  == 0 or i == len(strip)-1 :
			new_strip.append(strip[i])
		else:
			line = strip[i]
			for j in range(start,end):
				if  j != 0 or j != len(strip[i])-1:
					line = color_line(line, color,j)
			new_strip.append(line)
	return new_strip

def adjust_length(length, maxlength):
	percent = float(length)/float(maxlength)
	return int(percent * WINDOW_WIDTH)

def create_white_strip(length):
	strip = []
	strip.append(produce_line(length+2, black))
	for i in range(int(0.02*length)):
		line = produce_line(length+2,white)
		line = color_line(line,black,0)
		line = color_line(line, black,length+1)
		strip.append(line)
	strip.append(produce_line(length+2, black))
	return strip

def color_genome(location, chr_lengths, chr_names, genome, chr_colors):
	try:
		chromosome, start, end = location.split("_")
		chromosome = chr_names.index(chromosome)
		start = adjust_length(int(start),chr_lengths[chromosome])
		end = adjust_length(int(end), chr_lengths[chromosome])
		if end - start == 0:
			end = start+1
		genome[chromosome] = color_strip(genome[chromosome],chr_colors[chromosome],start,end)
	except:
		pass

chr_lengths = [
	248956422,
	242193529,
	198295559,
	190214555,
	181538259,
	170805979,
	159345973,
	145138636,
	138394717,
	133797422,
	135086622,
	133275309,
	114364328,
	107043718,
	101991189,
	90338345,
	83257441,
	80373285,
	58617616,
	64444167,
	46709983,
	50818468,
	156040895,
	57227415,
	16569]

chr_names = [
	"1",
	"2",
	"3",
	"4",
	"5",
	"6",
	"7",
	"8",
	"9",
	"10",
	"11",
	"12",
	"13",
	"14",
	"15",
	"16",
	"17",
	"18",
	"19",
	"20",
	"21",
	"22",
	"X",
	"Y",
	"MT"]

chr_colors = [
	[0, 150, 255],
	[0, 170, 235],
	[0, 190, 215],
	[0, 210, 195],
	[0, 230, 175],
	[0, 250, 155],
	[20, 250, 135],
	[40, 150, 115],
	[60, 150, 95],
	[80, 150, 75],
	[100, 150, 55],
	[120, 150, 35],
	[140, 150, 15],
	[160, 130, 15],
	[180, 110, 15],
	[210, 90, 15],
	[230, 70, 15],
	[250, 50, 15],
	[270, 30, 15],
	[150, 10, 255],
	[150, 150, 255],
	[150, 150, 255],
	[150, 150, 255],
	[150, 150, 255],
	[150, 150, 255]
]

genome = []
for i in range(len(chr_names)):
	genome.append(create_white_strip(WINDOW_WIDTH))

def display_locations(code, filename, chr_lengths, chr_names, genome, chr_colors, ref_number):
	with open(filename,"r") as f:
		for i, s in enumerate(f):
			s = s.rstrip().split("\t")
			if s[0:ref_number] == code.split("\t"):
				locations = s[ref_number:]
				for location in locations:
					print(location)
					color_genome(location, chr_lengths, chr_names, genome, chr_colors)
				break
	fix, ax = plt.subplots(len(genome),1)
	for i in range(len(genome)):
		ax[i].axis("off")
		ax[i].text(-40,15,chr_names[i])
		ax[i].imshow(genome[i])
	plt.suptitle(code)
	#plt.subplots_adjust(hspace = 2)
	plt.show()