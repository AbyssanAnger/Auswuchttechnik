from scipy.io import loadmat

dir = r"C:\Users\Sebastian\Desktop\Coding\Auswuchttechnik\awt_UE1_Aufgabe.mat"
data = loadmat(dir)
print(data.keys())
