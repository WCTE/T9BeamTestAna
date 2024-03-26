#!/usr/bin/python
# jk 24.4.2018

gParticles = {}

class MyParticle:
    def Set(self, M, z, name, tlatexname, texname):
        self._M = M
        self._z = z
        self._name = name
        self._tlatexname = tlatexname
        self._texname = texname

    def __init__(self, M, z, name, tlatexname, texname):
        self.Set(M, z, name, tlatexname, texname)
                
    def GetM(self): return self._M
    def GetZ(self): return self._z
    def GetName(self): return self._name
    def GetTLatexName(self): return self._tlatexname
    def GetTexName(self): return self._texname
    
# masses in MeV, chages in units of |e|
Particles = [ [ 938.2721, +1, 'Proton', 'p', 'p'],
              [ 0.511,    +1, 'Electron', 'e^{-}', 'e^{-}'],
              [ 0.511,    -1, 'Positron', 'e^{+}', 'e^{+}'],
              [ 105.6583716,  +1, 'Muon', '#mu', r'\mu'],
              [ 139.57,   +1, 'Pion', '#pi', r'\pi'],
              [ 493.677,   +1, 'Kaon', 'K^{#pm}', r'K^{\pm}'],
              [ 3727.,     +2, 'Alpha', '#alpha', r'\alpha'],
              [ 1875.6,   +1, 'Deuteron', 'D', 'D'],
              [ 3.01604928*931.494, +1, 'Tritium', 'T', 'T' ]
]

def MakeParticles():
    print('Making particles...')
    for particle in Particles:
        gParticles[particle[2]] = (MyParticle(particle[0], particle[1], particle[2], particle[3], particle[4]))

def PrintParticles():
    for part in gParticles:
        print(part)
