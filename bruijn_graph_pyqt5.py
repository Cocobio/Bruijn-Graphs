from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import networkx as nx
import random

import colorsys

import ast

import my_networkx as my_nx

class PrettyWidget(QWidget):

	NumButtons = ['updateGView']
	n = 5
	I = [0,1]
	K = [i for i in range(n) if i not in [0,1]]
	BasicConstrains = []
	AddedConstrains = []

	def __init__(self):


		super(PrettyWidget, self).__init__()		
		font = QFont()
		font.setPointSize(16)
		self.initUI()

	def initUI(self):

		self.setGeometry(100, 100, 800, 600)
		self.center()
		self.setWindowTitle('Bruijn Graph on Feedback Vertex Register')

		self.init_G()

		grid = QGridLayout()
		self.setLayout(grid)
		self.createVerticalArcsGroup()
		self.createVerticalGroupBox() 

		buttonLayout = QVBoxLayout()
		buttonLayout.addWidget(self.verticalGroupBox)

		self.figure = plt.figure()
		self.canvas = FigureCanvas(self.figure)	
		grid.addWidget(self.canvas, 0, 1, 9, 9)	
		# grid.addLayout(buttonLayout, 0, 0)

		self.settingDialog = QDialog(self)
		self.w = QWidget()
		self.w.setLayout(buttonLayout)
		# self.settingDialog.setCentralWidget(self.w)



		self.show()
		# self.settingDialog.show()
		self.w.show()

		self.updateGView()

	def createVerticalArcsGroup(self):
		self.arcListWidget = QListWidget()
		self.fillArcList()

	def createVerticalGroupBox(self):
		self.verticalGroupBox = QGroupBox()

		layout = QVBoxLayout()

		### modificacion de n
		layout_h = QHBoxLayout()
		label = QLabel()
		label.setText("n: ")
		layout_h.addWidget(label)

		lineEdit = QLineEdit()
		lineEdit.setObjectName("nLineEdit")
		lineEdit.setText(str(self.n))
		layout_h.addWidget(lineEdit)
		lineEdit.returnPressed.connect(self.changeN)

		layout.addLayout(layout_h)

		### Lista de nodos a ignorar y nodos a observar
		layout_h = QHBoxLayout()
		label = QLabel()
		label.setText("Ignore: ")
		layout_h.addWidget(label)

		lineEdit = QLineEdit()
		lineEdit.setObjectName("ILineEdit")
		lineEdit.setText(str(self.I))
		layout_h.addWidget(lineEdit)
		lineEdit.returnPressed.connect(self.modifyI)


		# layout_h = QHBoxLayout()
		label = QLabel()
		label.setText("Observe: ")
		layout_h.addWidget(label)

		lineEdit = QLineEdit()
		lineEdit.setObjectName("KLineEdit")
		lineEdit.setText(str(self.K))
		layout_h.addWidget(lineEdit)
		lineEdit.returnPressed.connect(self.modifyK)

		layout.addLayout(layout_h)

		button = QPushButton("updateGView")
		button.setObjectName("updateGView")
		layout.addWidget(button)
		layout.setSpacing(10)
		self.verticalGroupBox.setLayout(layout)
		button.clicked.connect(self.submitCommand)

		layout_h = QHBoxLayout()
		lineEdit = QLineEdit()
		lineEdit.setObjectName("fromLineEdit")
		layout_h.addWidget(lineEdit)
		lineEdit.returnPressed.connect(self.Add)

		label = QLabel()
		label.setText("->")
		layout_h.addWidget(label)
		
		lineEdit = QLineEdit()
		lineEdit.setObjectName("toLineEdit")
		layout_h.addWidget(lineEdit)
		lineEdit.returnPressed.connect(self.Add)

		for i in ["Add", "Remove"]:
			button = QPushButton(i)
			button.setObjectName(i)
			layout_h.addWidget(button)
			button.clicked.connect(self.submitCommand)

		layout.addLayout(layout_h)
		layout.addWidget(self.arcListWidget)

		self.arcListWidget.itemDoubleClicked.connect(self.deleteArc)

		button = QPushButton("List SCC")
		button.setObjectName("ListSCC")
		layout.addWidget(button)
		self.verticalGroupBox.setLayout(layout)
		button.clicked.connect(self.submitCommand)
		# self.verticalGroupBox.setLayout(layout)

		button = QPushButton("Run my function")
		button.setObjectName("runThisScript")
		layout.addWidget(button)
		self.verticalGroupBox.setLayout(layout)
		button.clicked.connect(self.submitCommand)

	def submitCommand(self):
		eval('self.' + str(self.sender().objectName()) + '()')


	def init_G(self):
		self.node_n = 2**self.n
		self.BasicConstrains = [(0,1),(self.node_n-1,self.node_n-2)]
		mask = (self.node_n-1)

		self.G = nx.DiGraph()
		self.node_labels = {}

		for i in range(self.node_n):
			self.G.add_node(i,subset=sum([(i>>j)&1 for j in range(self.n)]))

			arc = (i,(i<<1)&mask)
			self.G.add_edge(*arc)

			arc = (i,((i<<1)&mask)+1)
			self.G.add_edge(*arc)

			self.node_labels[i] = "{0:0{1}b}".format(i,self.n)

		self.I.sort()
		# Delete nodes that are out of bound
		while len(self.I)>0 and self.I[-1]>=self.n-1:
			self.I.pop()
		self.I_mask = sum([1<<i for i in self.I])
		self.n_mask = 1<<(self.n-1)
		cleaning_mask = self.n_mask|self.I_mask

		# Esto ayuda a poder eliminar cualquier combinacion de indices
		unique_x = set()
		x_to_color_id = {}

		for x in range(self.node_n):
			x_masked = x &~cleaning_mask

			if x_masked in unique_x:
				x_to_color_id[x] = x_to_color_id[x_masked]
			else:
				x_to_color_id[x] = len(unique_x)
				unique_x.add(x_masked)


		# print(arcs_labels)
		self.positive_rel = {}
		# negative_rel = {}

		for x in range(self.node_n):
			# negative_rel[x] = x^self.n_mask
			self.positive_rel[x&~self.I_mask] = [x] if (x&~self.I_mask) not in self.positive_rel else self.positive_rel[x&~self.I_mask] + [x]

		for arc in self.BasicConstrains:
			self.take_arc(*arc)

		for arc in self.AddedConstrains:
			self.take_arc(*arc)

		colors_n = self.node_n//2//(2**len(self.I))
		colors = []

		for i in range(colors_n):
			hsv_color = (i/colors_n,1,1)
			rgb_color = colorsys.hsv_to_rgb(*hsv_color)

			colors.append("#{0:0{3}x}{1:0{3}x}{2:0{3}x}".format(int(rgb_color[0]*255),int(rgb_color[1]*255),int(rgb_color[2]*255),2))

			# colors = ["#{0:0{3}x}{1:0{3}x}{2:0{3}x}".format(rgb_color[0],rgb_color[1],rgb_color[2],2) for i in range(colors_n)]
		self.node_colors = [colors[x_to_color_id[i & ~cleaning_mask]] for i in range(self.node_n)]

	def take_arc(self,x,y):
		adj = self.G.adj[x]

		if y not in adj:
			print("Error, contradiction on {0:0{2}b} -> {1:0{2}b}".format(x,y,self.n))
			return

		# Remove the other arc
		if len(adj)!=2:
			return
			
		self.G.remove_edge(x,y^1) ## Remove the other edge

		# Remove the negative rel arc
		if y in self.G.adj[x^self.n_mask]:
			self.take_arc(x^self.n_mask,y^1)

		mask = (self.node_n-1)
		for x_2 in self.positive_rel[x&~self.I_mask]:
			self.take_arc(x_2,((x_2<<1)|(y&1))&mask)

	def updateGView(self):
		self.figure.clf()
		############ Create 

		position = nx.multipartite_layout(self.G)

		############ Arcs label
		arcs_labels = {}
		for arc in self.G.edges():
			arcs_labels[arc] = arc[1]&1

		mandatory_arcs = []
		for v in range(self.node_n):
			if self.G.out_degree(v) == 1:
				mandatory_arcs.append((v,next(self.G.neighbors(v))))
				# G[v][next(G.neighbors(v))]['color'] = 'r'
			elif self.G.out_degree(v) == 0:
				raise("No outgoing arcs for: {0:0{1}b}".format(v,self.n))


		nx.draw_networkx_nodes(self.G, position, node_color=self.node_colors)#, labels=self.node_labels)
		nx.draw_networkx_labels(self.G, position, self.node_labels, font_size=8)

		##### Mandatory arcs
		nx.draw_networkx_edges(self.G, position, edgelist=mandatory_arcs, edge_color="tab:red", connectionstyle="arc3,rad=-0.21")

		##### Other arcs
		nx.draw_networkx_edges(self.G, position, edgelist=list(set([e for e in self.G.edges]) -set(mandatory_arcs)), connectionstyle="arc3,rad=-0.21")

		# nx.draw_networkx_edge_labels(self.G, position, edge_labels=arcs_labels)
		my_nx.my_draw_networkx_edge_labels(self.G, position, edge_labels=arcs_labels,rotate=False,rad = -0.21)
		self.canvas.draw_idle()

	def fillArcList(self):
		self.arcListWidget.clear()

		for i in range(len(self.BasicConstrains)):
			v,u = self.BasicConstrains[i]
			self.arcListWidget.insertItem(i, "({0:0{2}b},{1:0{2}b})".format(v,u,self.n))

		for i in range(len(self.AddedConstrains)):
			v,u = self.AddedConstrains[i]
			self.arcListWidget.insertItem(i+len(self.BasicConstrains),"({0:0{2}b},{1:0{2}b})".format(v,u,self.n))

	def center(self):
		qr = self.frameGeometry()
		cp = QDesktopWidget().availableGeometry().center()
		qr.moveCenter(cp)
		self.move(qr.topLeft())

	def Add(self):
		v = self.w.findChild(QLineEdit, "fromLineEdit").text()
		u = self.w.findChild(QLineEdit, "toLineEdit").text()

		if len(v)!=self.n or len(u)!=self.n:
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Critical)
			msg.setText("Error")
			msg.setInformativeText('Vertex information must be equal to the tags on the Graph.')
			msg.setWindowTitle("Error")
			msg.exec_()

			return

		v = int(v,2)
		u = int(u,2)

		if (v,u) in self.AddedConstrains:
			print("Arc already added ass Contrain")
			return
		if (v,u) in self.BasicConstrains:
			print("Arc is one of the basic constrains of Bruijn graphs")
			return
		if not self.G.has_edge(v,u):
			print("Arc does not exists on G")
			return
		if (self.G.has_edge(v,u) and self.G.out_degree(v)==1):
			print("Arc is already the only posible outgoing arc for {0:0{1}b}".format(v,self.n))
			return

		self.AddedConstrains.append((v,u))
		self.take_arc(v,u)
		self.fillArcList()
		self.updateGView()


	def Remove(self):
		v = self.w.findChild(QLineEdit, "fromLineEdit").text()
		u = self.w.findChild(QLineEdit, "toLineEdit").text()

		if len(v)!=self.n or len(u)!=self.n:
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Critical)
			msg.setText("Error")
			msg.setInformativeText('Vertex information must be equal to the tags on the Graph.')
			msg.setWindowTitle("Error")
			msg.exec_()

			return

		v = int(v,2)
		u = int(u,2)

		if (v,u) not in self.AddedConstrains:
			print("Constrain not on the list.")
			return

		self.AddedConstrains.remove((v,u))

		self.init_G()
		self.fillArcList()
		self.updateGView()

	def ListSCC(self):
		print("Printing the strongly connected components on the Bruijn Graph!")
		i = 1
		for component in nx.strongly_connected_components(self.G):
			print("Component ",i)

			# print("\tFirst node of component {0}: {1:0{2}b}".format(i,next(iter(component)), self.n))
			for v in component:
				print("\t{0:0{1}b}".format(v,self.n))

			i += 1

	def changeN(self):
		try:
			new_n = int(self.w.findChild(QLineEdit, "nLineEdit").text())
			self.n = new_n
		except:
			print("Can't parse text as integer")
			return

		self.init_G()
		# self.fillArcList()
		self.updateGView()

		self.w.findChild(QLineEdit, "ILineEdit").setText(str(self.I))

	def modifyI(self):
		data = self.w.findChild(QLineEdit, "ILineEdit").text()

		old_I = [i for i in self.I]

		try:
			self.I = ast.literal_eval(data)
			self.I = list(set(self.I))
		except:
			print("Can't parse text as Python list")
			return

		if type(self.I) != type(old_I):
			self.I = old_I
			print("No list found on line edit.")
			return

		self.K = [i for i in range(self.n) if i not in self.I]
		self.init_G()
		# self.fillArcList()
		self.updateGView()

		self.w.findChild(QLineEdit, "ILineEdit").setText(str(self.I))
		self.w.findChild(QLineEdit, "KLineEdit").setText(str(self.K))

	def modifyK(self):
		data = self.w.findChild(QLineEdit, "KLineEdit").text()

		old_K = [i for i in self.K]

		try:
			self.K = ast.literal_eval(data)
			self.K = list(set(self.K+[self.n-1]))
		except:
			print("Can't parse text as Python list")
			return

		if type(self.K) != type(old_K):
			self.K = old_K
			print("No list found on line edit.")
			return

		self.K.sort()
		while len(self.K)>0 and self.K[-1] >= self.n:
			self.K.pop()

		self.I = [i for i in range(self.n) if i not in self.K]

		self.init_G()
		self.updateGView()

		self.w.findChild(QLineEdit, "ILineEdit").setText(str(self.I))
		self.w.findChild(QLineEdit, "KLineEdit").setText(str(self.K))

	def deleteArc(self,item):
		arc = item.text()[1:-1].split(',')
		print(arc)
		arc = (int(arc[0],2), int(arc[1],2))

		if arc in self.BasicConstrains:
			return

		self.AddedConstrains.remove(arc)

		self.init_G()
		self.fillArcList()
		self.updateGView()

	def runThisScript(self):
		## The graph is saved in self.G
		## To force an arc, use the function self.take_arc
		## Any questions, ask Ignacio :)
		print("This part is running!")


if __name__ == '__main__':

	import sys  
	app = QApplication(sys.argv)
	app.aboutToQuit.connect(app.deleteLater)
	app.setStyle(QStyleFactory.create("gtk"))
	screen = PrettyWidget() 
	screen.show()   
	sys.exit(app.exec_())