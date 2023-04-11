import string
from molgrafik import *
from atomdict import *
"""
<formel>::= <mol> \n
<mol>   ::= <group> | <group><mol>
<group> ::= <atom> |<atom><num> | (<mol>) <num>
<atom>  ::= <LETTER> | <LETTER><letter>
<LETTER>::= A | B | C | ... | Z
<letter>::= a | b | c | ... | z
<num>   ::= 2 | 3 | 4 | ...
"""

atomlist = [
  "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si",
  "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co",
  "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
  "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",
  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
  "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au",
  "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",
  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db",
  "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Fl", "Lv"
]


class Syntaxfel(Exception):
  pass


class Ruta:

  def __init__(self, atom="( )", num=1):
    self.atom = atom
    self.num = num
    self.next = None
    self.down = None

  def __str__(self):
    return " Vikt atom {}*{}".format(self.atom,self.num)


def readmolekyl(q):
  mol = readgroup(q)
  peek = q.peek()
  if peek == ".":
    return mol
  elif peek == ")":
    return mol
  mol.next = readmolekyl(q)


def readgroup(q):
  rutan = Ruta()
  peek = q.peek()
  if peek == "(":
    q.dequeue()
    rutan.down = readmolekyl(q)
    peek = q.peek()
    if peek == ")":
      q.dequeue()
      peek = q.peek()
      if not peek.isdigit():
        raise Syntaxfel("Saknad siffra vid radslutet ")
      rutan.num = readnum(q)
      peek = q.peek()
      if peek.isalpha() or peek == "(":
        rutan.next = readmolekyl(q)
      return rutan
    raise Syntaxfel("Saknad högerparentes vid radslutet ")
  elif peek.isalpha():
    rutan.atom = readatom(q)
    peek = q.peek()
    if peek.isdigit():
      rutan.num = readnum(q)
    peek = q.peek()
    if peek.isalpha() or peek == "(":
      rutan.next = readmolekyl(q)
    return rutan
  elif peek == ".":
    return rutan
  else:
    raise Syntaxfel("Felaktig gruppstart vid radslutet ")


def readatom(q):
  l1 = readLETTER(q)
  peek = q.peek()
  if not peek.islower():
    if l1 in atomlist:
      return l1
    raise Syntaxfel("Okänd atom vid radslutet ")
  l2 = readletter(q)
  if l1 + l2 in atomlist:
    return l1 + l2
  raise Syntaxfel("Okänd atom vid radslutet ")


def readLETTER(q):
  L = q.dequeue()
  string_upper = string.ascii_uppercase
  LETTER = list(string_upper)
  if L in LETTER:
    return L
  raise Syntaxfel("Saknad stor bokstav vid radslutet " + L)


def readletter(q):
  l = q.dequeue()
  string_lower = string.ascii_lowercase
  letter = list(string_lower)
  if l in letter:
    return l


def readnum(q):
  n = q.dequeue()
  num = ""
  num += n
  if int(n) >= 2:
    peek = q.peek()
    while peek.isdigit():
      n = q.dequeue()
      num += n
      peek = q.peek()
    return int(num)
  elif int(n) == 1:
    peek = q.peek()
    if not peek.isdigit():
      raise Syntaxfel("För litet tal vid radslutet ")
    peek = q.peek()
    while peek.isdigit():
      n = q.dequeue()
      num += n
      peek = q.peek()
    return int(num)
  else:
    raise Syntaxfel("För litet tal vid radslutet ")


def weight(mol):
  weightdict = skapaAtomdict()
  molweight = 0 # tror vi måsta ha såhär här annars säger den att den inte är definierad när den returneras
  if mol.atom in weightdict:
    molweight = weightdict.get(mol.atom)
    molweight *= mol.num 
    print(mol, "=", molweight)
  
  if not mol.down is None:
    # print("down")
    total = 0 # ny total för det i paranteserna
    total += weight(mol.down) #lägg på det här som också kör alla next o lägger på o returnerar
    total *= mol.num # gångra med parantes num
    print(" I parentes *", mol.num, "=", total)
    molweight += total # lägg på på molweight
  if not mol.next is None:
    # print("next")
    molweight += weight(mol.next)
  return molweight


def printQueue(q):
  string = ""
  while not q.isEmpty():
    d = q.dequeue()
    if d != ".":
      string += d
  return string


def readformel(q):
  mol = readmolekyl(q)
  if q.peek() == ")":
    raise Syntaxfel("Felaktig gruppstart vid radslutet ")
  return mol


def main():
  while True:
    mg = Molgrafik()
    molekyl = input()
    if molekyl == "#":
      break
    q = LinkedQ()
    m = list(molekyl)
    for delar in m:
      q.enqueue(delar)
    q.enqueue(".")
    try:
      mol = readformel(q)
      mg.show(mol)     #HÄR ÄR ANROPET
      print("Formeln är syntaktiskt korrekt")
      print("Totala vikten:", weight(mol))
    except Syntaxfel as fel:
      print(str(fel) + printQueue(q))


#_______________________________________________________________
#LinkedQ


class Node:

  def __init__(self, value, next=None):
    self.value = value
    self.next = next


class LinkedQ:

  def __init__(self):
    self._first = None
    self._last = None

  def enqueue(self, x):
    new = Node(x)
    if self.isEmpty():
      self._first = new
      self._last = self._first

    else:
      self._last.next = new
      self._last = self._last.next

  def dequeue(self):
    if self.isEmpty():
      return None
    else:
      x = self._first.value
      self._first = self._first.next
      return x

  def isEmpty(self):
    """kollar om kon är tom, behovs till enqueue och dequeue"""
    if self._first is None:
      return True
    else:
      return False

  def peek(self):
    return self._first.value


if __name__ == '__main__':
  main()
