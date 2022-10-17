from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import networkx as nx
import random

import colorsys

import my_networkx as my_nx

class PrettyWidget(QWidget):

	NumButtons = ['plot1','plot2', 'plot3']
	n = 5
	I = [0,1]

	def __init__(self):


		super(PrettyWidget, self).__init__()		
		font = QFont()
		font.setPointSize(16)
		self.initUI()

	def initUI(self):

		self.setGeometry(100, 100, 800, 600)
		self.center()
		self.setWindowTitle('Bruijn Graph on Feedback Vertex Register')

		grid = QGridLayout()
		self.setLayout(grid)
		self.createVerticalGroupBox() 

		buttonLayout = QVBoxLayout()
		buttonLayout.addWidget(self.verticalGroupBox)

		self.figure = plt.figure()
		self.canvas = FigureCanvas(self.figure)	
		grid.addWidget(self.canvas, 0, 1, 9, 9)	
		grid.addLayout(buttonLayout, 0, 0)

		self.show()


	def createVerticalGroupBox(self):
		self.verticalGroupBox = QGroupBox()

		layout = QVBoxLayout()
		for i in  self.NumButtons:
			button = QPushButton(i)
			button.setObjectName(i)
			layout.addWidget(button)
			layout.setSpacing(10)
			self.verticalGroupBox.setLayout(layout)
			button.clicked.connect(self.submitCommand)

	def submitCommand(self):
		eval('self.' + str(self.sender().objectName()) + '()')



	def plot1(self):
		self.figure.clf()
		ax1 = self.figure.add_subplot(211)
		x1 = [i for i in range(100)]
		y1 = [i**0.5 for i in x1]
		ax1.plot(x1, y1, 'b.-')

		ax2 = self.figure.add_subplot(212)
		x2 = [i for i in range(100)]
		y2 = [i for i in x2]
		ax2.plot(x2, y2, 'b.-')
		self.canvas.draw_idle()

	def plot2(self):
		self.figure.clf()
		ax3 = self.figure.add_subplot(111)
		x = [i for i in range(100)]
		y = [i**0.5 for i in x]
		ax3.plot(x, y, 'r.-')
		ax3.set_title('Square Root Plot')
		self.canvas.draw_idle()

	def plot3(self):
		self.figure.clf()
		node_n = 2**self.n
		mask = (node_n-1)

		G = nx.DiGraph()
		node_labels = {}

		for i in range(node_n):
			G.add_node(i,subset=sum([(i>>j)&1 for j in range(self.n)]))

			arc = (i,(i<<1)&mask)
			G.add_edge(*arc)

			arc = (i,((i<<1)&mask)+1)
			G.add_edge(*arc)

			node_labels[i] = "{0:0{1}b}".format(i,self.n)

		self.I.sort()
		# Delete nodes that are out of bound
		while len(self.I)>0 and self.I[-1]>=self.n-1:
			self.I.pop()
		I_mask = sum([1<<i for i in self.I])
		n_mask = 1<<(self.n-1)
		cleaning_mask = n_mask|I_mask

		# Esto ayuda a poder eliminar cualquier combinacion de indices
		unique_x = set()
		x_to_color_id = {}

		for x in range(node_n):
			x_masked = x &~cleaning_mask

			if x_masked in unique_x:
				x_to_color_id[x] = x_to_color_id[x_masked]
			else:
				x_to_color_id[x] = len(unique_x)
				unique_x.add(x_masked)

		colors_n = node_n//2//(2**len(self.I))

		colors = []

		for i in range(colors_n):
			hsv_color = (i/colors_n,1,1)
			rgb_color = colorsys.hsv_to_rgb(*hsv_color)

			colors.append("#{0:0{3}x}{1:0{3}x}{2:0{3}x}".format(int(rgb_color[0]*255),int(rgb_color[1]*255),int(rgb_color[2]*255),2))

			# colors = ["#{0:0{3}x}{1:0{3}x}{2:0{3}x}".format(rgb_color[0],rgb_color[1],rgb_color[2],2) for i in range(colors_n)]
		node_colors = [colors[x_to_color_id[i & ~cleaning_mask]] for i in range(node_n)]


		# print(arcs_labels)
		positive_rel = {}
		negative_rel = {}

		for x in range(node_n):
			negative_rel[x] = x^n_mask
			positive_rel[x&~I_mask] = [x] if (x&~I_mask) not in positive_rel else positive_rel[x&~I_mask] + [x]

		def take_arc(x,y):
			adj = G.adj[x]

			if y not in adj:
				print("Error, contradiction on {0:0{2}b} -> {1:0{2}b}".format(x,y,self.n))
				return

			# Remove the other arc
			if len(adj)!=2:
				return
			
			G.remove_edge(x,y^1) ## Remove the other edge

			# Remove the negative rel arc
			if y in G.adj[x^n_mask]:
				take_arc(x^n_mask,y^1)

			for x_2 in positive_rel[x&~I_mask]:
				take_arc(x_2,((x_2<<1)|(y&1))&mask)

		############ Propagate constrain 0 -> 1
		take_arc(0,1)

		############ Propagate constrain 1^n -> 1^(n-1)0		
		take_arc(node_n-1,node_n-2)

		# take_arc(7,14)

		############ Create 
		position = nx.multipartite_layout(G)

		############ Arcs label
		arcs_labels = {}
		for arc in G.edges():
			arcs_labels[arc] = arc[1]&1

		nx.draw_networkx_nodes(G, position, node_color=node_colors)#, labels=node_labels)
		nx.draw_networkx_labels(G, position, node_labels, font_size=8)
		nx.draw_networkx_edges(G, position, connectionstyle="arc3,rad=-0.21")
		# nx.draw_networkx_edge_labels(G, position, edge_labels=arcs_labels)
		my_nx.my_draw_networkx_edge_labels(G, position, edge_labels=arcs_labels,rotate=False,rad = -0.21)
		self.canvas.draw_idle()

	def center(self):
		qr = self.frameGeometry()
		cp = QDesktopWidget().availableGeometry().center()
		qr.moveCenter(cp)
		self.move(qr.topLeft())

if __name__ == '__main__':

	import sys  
	app = QApplication(sys.argv)
	app.aboutToQuit.connect(app.deleteLater)
	app.setStyle(QStyleFactory.create("gtk"))
	screen = PrettyWidget() 
	screen.show()   
	sys.exit(app.exec_())