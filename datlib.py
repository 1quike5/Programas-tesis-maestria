import numpy as np
import math 
import os.path

class muestras:
    def __init__(self, ruta, nombre, n):
        self.rutaL = []
        self.rutaT = []
        for i in range(n):
            self.rutaL.append(ruta + "/" + nombre + str(i+1) + "/luminiscencia.txt")
            self.rutaT.append(ruta + "/" + nombre + str(i+1) + "/transmitancia.txt")
            self.rutaP = ruta + "/patron.txt"
            self.setDataSize()
            self.dataLoad()
            self.graphData()
    def setDataSize(self):
        n = np.zeros(len(self.rutaL))
        dataSize = np.zeros(len(self.rutaL))
        for i in range(len(n)):
            j = 0
            k = 0
            with open(self.rutaL[i], 'r') as f:
                for line in f:
                    j = j + 1
                    if (line.find(".Foto") or line.find("-- Fecha")) >= 0:
                        dataSize[i] = j - k - 2
                        k = j
            n[i] = int(j/(dataSize[i]+2))
        self.n = n
        self.dataSize = dataSize
    def dataLoad(self):
        self.wavelength = np.loadtxt(self.rutaL[0], skiprows = 2, max_rows = int(self.dataSize[0]), delimiter = '\t', usecols = 0)
        self.patron = np.loadtxt(self.rutaP, skiprows = 2, max_rows = int(self.dataSize[0]), delimiter = '\t', usecols = 1)
        datosG = np.empty([len(self.n), int(min(self.n)), int(self.dataSize[0])])
        datosGT = np.empty([len(self.n), int(min(self.n)), int(self.dataSize[0])])
        datosGA = np.empty([len(self.n), int(min(self.n)), int(self.dataSize[0])])
        for j in range(len(self.n)):
            datos = np.array([])
            datosT = np.array([])
            
            for i in range(int(self.n[j])):
                lectura = np.loadtxt(self.rutaL[j], skiprows = 2 + 2*i + int(self.dataSize[j])*i, max_rows = int(self.dataSize[j]), delimiter = '\t', usecols = 1)
                lecturaT = np.loadtxt(self.rutaT[j], skiprows = 2 + 2*i + int(self.dataSize[j])*i, max_rows = int(self.dataSize[j]), delimiter = '\t', usecols = 1)
                lecturaT = lecturaT/self.patron
                for k in range(lecturaT.size):
                    if math.isnan(lecturaT[int(k)]):
                        lecturaT[k] = 0.00000001
                datos = np.concatenate([datos, lectura])
                datosT = np.concatenate([datosT, lecturaT])
                
            absorbancia = - np.log(datosT)
            absorbancia = absorbancia.reshape(int(self.n[j]), int(self.dataSize[j]))

            datos = datos.reshape(int(self.n[j]), int(self.dataSize[j]))
            datosT = datosT.reshape(int(self.n[j]), int(self.dataSize[j]))
            datosG[j,:,:] = datos[0:int(min(self.n)), :]
            datosGT[j,:,:] = datosT[0:int(min(self.n)), :]
            datosGA[j,:,:] = absorbancia[0:int(min(self.n)), :]
        size = len(self.wavelength)
        cut = 0
        upper_cut = size
        for i in range(size):
            if self.wavelength[i] >= 425:
                cut = i
                break
        for i in range(size):
            if self.wavelength[i] >= 684:
                upper_cut = i
                break
        self.wavelength = self.wavelength[cut:upper_cut]
        datosG = datosG[:,:,cut:upper_cut]
        datosGT = datosGT[:,:,cut:upper_cut]
        datosGA = datosGA[:,:,cut:upper_cut]
        self.datos = datosG
        self.datosT = datosGT
        self.absorbancia = datosGA
    def graphData(self):
        self.offset = self.datos[:,0,0]
        maximo = []
        for i in range(len(self.n)):
            maximo.append(max(self.datos[i,0,:]))
        self.maximo = max(maximo)
        self.maximos = maximo
        self.factor = self.maximo/maximo
        self.pfactor = np.ones(len(self.n))        
        
class transmitance:
    def __init__(self, ruta, nombre, modo, normalize):
        self.modo = modo
        self.normalize = normalize
        self.M = 2
        if modo == "comparar":
            self.rutaC = ruta + "/" + nombre + "C.txt"
            self.rutapC = ruta + "/patronC.txt"
            self.M = 3
        self.rutaN = ruta + "/" + nombre + "1.txt"
        self.rutaO = ruta + "/" + nombre + "2.txt"
        self.rutapN = ruta + "/patron1.txt"
        self.rutapO = ruta + "/patron2.txt"
        self.setDataSize()
        self.dataLoad()
        self.graphData()
    
    def setDataSize(self):
        n = 1
        dataSize = 0
        j = 0
        k = 0
        with open(self.rutaN, 'r') as f:
            for line in f:
                j = j + 1
                if (line.find(".Foto") or line.find("-- Fecha")) >= 0:
                    dataSize = j - k - 2
                    k = j
        n = int(j/(dataSize+2))
        self.n = n
        self.dataSize = dataSize
        if self.modo == "comparar":
            n = 1
            dataSize = 0
            j = 0
            k = 0
            with open(self.rutaC, 'r') as f:
                for line in f:
                    j = j + 1
                    if (line.find(".Foto") or line.find("-- Fecha")) >= 0:
                        dataSize = j - k - 2
                        k = j
            n = int(j/(dataSize+2))
            self.n = n
            self.dataSize = dataSize
        
    def dataLoad(self):
        self.wavelength = np.loadtxt(self.rutaN, skiprows = 2, max_rows = int(self.dataSize), delimiter = '\t', usecols = 0)
        self.patronN = np.loadtxt(self.rutapN, skiprows = 2, max_rows = int(self.dataSize), delimiter = '\t', usecols = 1)
        self.patronO = np.loadtxt(self.rutapO, skiprows = 2, max_rows = int(self.dataSize), delimiter = '\t', usecols = 1)
        
        if self.modo == "comparar":
            self.patronC = np.loadtxt(self.rutapC, skiprows = 2, max_rows = int(self.dataSize), delimiter = '\t', usecols = 1)
            self.datosC = np.array([])
        datosG = np.empty([self.M, self.n, int(self.dataSize)])
        datosA = np.empty([self.M, self.n, int(self.dataSize)])
        datosN = np.array([])
        datosO = np.array([])
        
        
        for i in range(int(self.n)):
            lecturaN = np.loadtxt(self.rutaN, skiprows = 2 + 2*i + int(self.dataSize)*i, max_rows = int(self.dataSize), delimiter = '\t', usecols = 1)
            lecturaN = lecturaN/self.patronN
            
            lecturaO = np.loadtxt(self.rutaO, skiprows = 2 + 2*i + int(self.dataSize)*i, max_rows = int(self.dataSize), delimiter = '\t', usecols = 1)
            lecturaO = lecturaO/self.patronO
            
            for k in range(lecturaN.size):
                if math.isnan(lecturaN[int(k)]):
                    lecturaN[k] = 0.00000001
            for k in range(lecturaO.size):
                if math.isnan(lecturaO[int(k)]):
                    lecturaO[k] = 0.00000001
                
            datosN = np.concatenate([datosN, lecturaN])
            datosO = np.concatenate([datosO, lecturaO])
            
            if self.modo == "comparar":
                self.lecturaC = np.loadtxt(self.rutaC, skiprows = 2 + 2*i + int(self.dataSize)*i, max_rows = int(self.dataSize), delimiter = '\t', usecols = 1)
                self.lecturaC = self.lecturaC/self.patronC
                
                for k in range(self.lecturaC.size):
                    if math.isnan(self.lecturaC[int(k)]):
                        self.lecturaC[k] = 0.00000001

                self.datosC = np.concatenate([self.datosC, self.lecturaC])

        
        
        absorbanciaN = 2.0 - np.log10(datosN*100)
        absorbanciaO = 2.0 - np.log10(datosO*100)
        absorbanciaN = absorbanciaN.reshape(int(self.n), int(self.dataSize))
        absorbanciaO = absorbanciaO.reshape(int(self.n), int(self.dataSize))
        datosA[0,:,:] = absorbanciaN[0:self.n, :]
        datosA[1,:,:] = absorbanciaO[0:self.n, :]
        
        
        datosN = datosN.reshape(int(self.n), int(self.dataSize))
        datosO = datosO.reshape(int(self.n), int(self.dataSize))
        datosG[0,:,:] = datosN[0:self.n, :]
        datosG[1,:,:] = datosO[0:self.n, :]
        
        if self.modo == "comparar":
            self.absorbanciaC = 2.0 - np.log10(self.datosC*100)
            datosA[2,:,:] = self.absorbanciaC[0:self.n, :]
            
            self.datosC = self.datosC.reshape(int(self.n), int(self.dataSize))
            datosG[2,:,:] = self.datosC[0:self.n, :]
            
            
            
        size = len(self.wavelength)
        cut = 0
        upper_cut = size
        for i in range(size):
            if self.wavelength[i] >= 425:
                cut = i
                break
        for i in range(size):
            if self.wavelength[i] >= 684:
                upper_cut = i
                break
        self.wavelength = self.wavelength[cut:upper_cut]
        
        datosG = datosG[:,:,cut:upper_cut]
        datosA = datosA[:,:,cut:upper_cut]
        if self.normalize:
            for i in range(self.M):
                norm = max(datosG[i,0,:]) 
                for j in range(self.n):
                    
                    datosG[i,j,:] = datosG[i,j,:]/norm
                    datosA[i,j,:] = - np.log10(datosG[i,j,:])
        
        
        self.datos = datosG
        self.absorbancia = datosA
        
    def graphData(self):
        self.offset = self.datos[:,0,0]
        maximo = []
        for i in range(self.M):
            maximo.append(max(self.datos[i,0,:]))
        self.maximo = max(maximo)
        self.maximos = maximo
        self.factor = self.maximo/maximo
        self.pfactor = np.ones(self.M)   
        
        self.offsetA = self.absorbancia[:,0,0]
        maximoA = []
        for i in range(self.M):
            maximoA.append(max(self.absorbancia[i,0,:]))
        self.maximoA = max(maximoA)
        self.maximosA = maximoA
        self.factorA = self.maximoA/maximoA
        self.pfactorA = np.ones(self.M)   
        
class comparison:
    def __init__(self, ruta, nombre):
        self.ruta1 = ruta + "/" + nombre + "1.txt"
        self.ruta2 = ruta + "/" + nombre + "2.txt"
        
        if not os.path.exists(self.ruta1):
            self.ruta1 = ruta + "/" + nombre + ".txt"
            self.ruta2 = ruta + "/" + nombre + ".txt"
        
        self.M = 2
        self.setDataSize()
        self.dataLoad()
        
    def setDataSize(self):
        n = 1
        dataSize = 0
        j = 0
        k = 0
        findCount = 0
        with open(self.ruta1, 'r') as f:
            for line in f:
                j = j + 1
                if (line.find(".Foto") or line.find("-- Fecha")) >= 0:
                    dataSize = j - k - 2
                    k = j
                    findCount += 1
                    
        if findCount <= 1 :
            dataSize = j - k - 2
        n = int(j/(dataSize+2))
        self.n = n
        self.dataSize = dataSize
    
    def dataLoad(self):
        self.wavelength = np.loadtxt(self.ruta1, skiprows = 2, max_rows = int(self.dataSize), delimiter = '\t', usecols = 0)
        
        datosG = np.empty([self.M, self.n, int(self.dataSize)])
        dif_espectro = np.empty([self.n, int(self.dataSize)])
        
        datos1 = np.array([])
        datos2 = np.array([])
        
        for i in range(int(self.n)):
            lectura1 = np.loadtxt(self.ruta1, skiprows = 2 + 2*i + int(self.dataSize)*i, max_rows = int(self.dataSize), delimiter = '\t', usecols = 1)
            lectura2 = np.loadtxt(self.ruta2, skiprows = 2 + 2*i + int(self.dataSize)*i, max_rows = int(self.dataSize), delimiter = '\t', usecols = 1)
            
            for k in range(lectura1.size):
                if math.isnan(lectura1[int(k)]):
                    lectura1[k] = 0.00000001
            for k in range(lectura2.size):
                if math.isnan(lectura2[int(k)]):
                    lectura2[k] = 0.00000001
                
            datos1 = np.concatenate([datos1, lectura1])
            datos2 = np.concatenate([datos2, lectura2])
        
        
        datos1 = datos1.reshape(int(self.n), int(self.dataSize))
        datos2 = datos2.reshape(int(self.n), int(self.dataSize))
        datosG[0,:,:] = datos1[0:self.n, :]
        datosG[1,:,:] = datos2[0:self.n, :]
        for i in range(self.n):
            dif_espectro[i,:] = datosG[0,i,:]/datosG[1,i,:] 
            
            
            
        size = len(self.wavelength)
        cut = 0
        upper_cut = size
        for i in range(size):
            if self.wavelength[i] >= 425:
                cut = i
                break
        for i in range(size):
            if self.wavelength[i] >= 684:
                upper_cut = i
                break
        self.wavelength = self.wavelength[cut:upper_cut]
        
        datosG = datosG[:,:,cut:upper_cut]
        dif_espectro = dif_espectro[:, cut:upper_cut]
        maxi = np.ones(self.M)
        for i in range(self.M):
            maxi[i] = max(datosG[i,0,:])
            
        for i in range(self.M):
            for j in range(self.n):
                datosG[i,j,:] = datosG[i,j,:]/maxi[i]
        for i in range(self.n):
            dif_espectro[i,:] = datosG[0,i,:]/datosG[1,i,:] 
        self.espectro = dif_espectro
        self.datos = datosG
        
class absorption:
    def __init__(self, ruta, carpeta, numero, nombre):
        self.ruta = ruta + "/" + carpeta + numero + "/" + nombre + ".txt"
        self.rutaP = ruta + "/patron.txt"
        self.setDataSize()
        self.dataLoad()
        self.graphData()
    
    def setDataSize(self):
        n = 1
        dataSize = 0
        j = 0
        k = 0
        with open(self.ruta, 'r') as f:
            for line in f:
                j = j + 1
                if (line.find(".Foto") or line.find("-- Fecha")) >= 0:
                    dataSize = j - k - 2
                    k = j
        n = int(j/(dataSize+2))
        self.n = n
        self.dataSize = dataSize
        
    def dataLoad(self):
        self.wavelength = np.loadtxt(self.ruta, skiprows = 2, max_rows = int(self.dataSize), delimiter = '\t', usecols = 0)
        self.patron = np.loadtxt(self.rutaP, skiprows = 2, max_rows = int(self.dataSize), delimiter = '\t', usecols = 1)
        
        datosA = np.empty([self.n, int(self.dataSize)])
        datos = np.array([])
        
        
        for i in range(int(self.n)):
            lectura = np.loadtxt(self.ruta, skiprows = 2 + 2*i + int(self.dataSize)*i, max_rows = int(self.dataSize), delimiter = '\t', usecols = 1)
            lectura = lectura/self.patron
            
            for k in range(lectura.size):
                if math.isnan(lectura[int(k)]):
                    lectura[k] = 0.00000001
                
            datos = np.concatenate([datos, lectura])

        
        
        absorbancia = 2.0 - np.log10(datos*100)
        absorbancia = absorbancia.reshape(int(self.n), int(self.dataSize))
        datosA = absorbancia[0:self.n, :]
        
        
        datos = datos.reshape(int(self.n), int(self.dataSize))
        
            
        size = len(self.wavelength)
        cut = 0
        upper_cut = size
        for i in range(size):
            if self.wavelength[i] >= 425:
                cut = i
                break
        for i in range(size):
            if self.wavelength[i] >= 684:
                upper_cut = i
                break
        cut = 0
        self.wavelength = self.wavelength[cut:upper_cut]
        
        datosA = datosA[:,cut:upper_cut]
        
        self.absorbancia = datosA
        self.datos = datos[:,cut:upper_cut]
        
    def graphData(self):
        self.offset = self.datos[0,0]
        maximo = []
        maximo.append(max(self.datos[0,:]))
        self.maximo = max(maximo)
        self.maximos = maximo
        self.factor = self.maximo/maximo
        
        self.offsetA = self.absorbancia[0,0]
        maximoA = []
        
        maximoA.append(max(self.absorbancia[0,:]))
        self.maximoA = max(maximoA)
        self.maximosA = maximoA
        self.factorA = self.maximoA/maximoA
def lorentziana(x, w, A, x0):
    f = (2.0*A/np.pi)*(w/(4.0*(x-x0)**2.0 + w*w))
    return f
def lorentziana_mult(x, y0, *args):
    peaks = int(len(args)/3)
    w  = args[0:peaks]
    A  = args[peaks:2*peaks]
    x0 = args[2*peaks:3*peaks]
    f = np.zeros(len(x))
    for i in range(peaks):
        f += (2.0*A[i]/np.pi)*(w[i]/(4.0*(x-x0[i])**2.0 + w[i]**2.0))
    return f + y0
def plotsLorentz(x, n, *args):
    peaks = int((len(args) - 1)/3)

    w  = args[1 + n]
    A  = args[1 + peaks + n]
    x0 = args[1 + 2*peaks + n]

    f = (2.0*A/np.pi)*(w/(4.0*(x-x0)**2.0 + w**2.0))
    return f
def cut_data(x, y, bottom, up):
    index = 0
    index2 = 0 
    found = False
    for i in range(len(x)):
        if x[i] > bottom and not found:
            index = i
            found = True
        if x[i] > up:
            index2 = i
            break
    if index2 == 0:
        index2 = len(x)-1
    return y[:,index:index2], index, index2