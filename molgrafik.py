# Henrik Eriksson 2003-11-17
# revision: Magnus Rosell 2006-08-10
# revision: Joel W m fl 2010-10-27
# revision: Linda 2011-09-27
# revision: Linda 2012-10-03 (Tkinter/tkinter)

def version():
    """Kolla vilken version av Python vi har"""
    import sys
    datatyp = type(sys.version_info)
    if datatyp == type(()):
        version = sys.version_info[0]
    else:
        version = sys.version_info.major
    return version

if version() == 3:
    from tkinter import *
else:
    from Tkinter import *


class Ruta:
    def __init__(self, atom = "()", num = 1):
        self.atom = atom
        self.num = num
        self.next=None
        self.down=None


class Molgrafik:

    def __init__(self):
        self.root=None
        self.stor=("Courier",18,"bold")
        self.liten=("Courier",14,"bold")

    def ram(self,master,sidan):
        """Returnerar en ram. Parametrar: master (grafikfonstret), sidan (vilken sida den ska ligga mot, t ex LEFT) """
        ramen=Frame(master,bg="white")
        ramen.pack(side=sidan,fill=BOTH)
        return ramen

    def atomruta(self,master,namn,num):
        """Ritar en atomruta. Parametrar: master (grafikfonstret), namn (atomnamnet), num (antal atomer) """
        ruta=Frame(master,bg="yellow",borderwidth=2,relief=GROOVE)
        ruta.pack(side=LEFT)
        atom=Frame(ruta,bg="yellow")
        atom.pack(side=LEFT)
        Label(atom,text=namn,font=self.stor,bg="yellow").pack()
        Frame(atom, height=5, bg="yellow").pack()
        if num>1:
            Label(ruta, text=str(num),font=self.liten,bg="yellow").pack(side=BOTTOM)

    def streck(self,master):
        """ Ritar ett streck. Parametrar: master (grafikfonstret) """
        strecket=Frame(master)
        strecket.pack(side=LEFT,fill=BOTH,expand=True)
        Frame(strecket,bg="white",height=20).pack(fill=X)
        Frame(strecket,bg="red",height=4,width=25).pack(fill=X)
        Frame(strecket,bg="white").pack(fill=BOTH,expand=1)

    def stolpe(self,master):
        """ Ritar en stolpe. Parametrar: master (grafikfonstret) """
        hela=self.ram(master,TOP)
        stolpen=self.ram(hela,LEFT)
        Frame(stolpen,bg="white",width=15).pack(side=LEFT)
        Frame(stolpen,bg="red",width=4,height=25).pack(side=LEFT)
        Frame(hela,bg="white").pack(fill=BOTH,expand=1)

    def picture(self,master,p):
        """ Ritar bilden. Parametrar: master (grafikfonstret), p (referens till datstrukturen som ska ritas) """
        if p is None: return
        storruta=self.ram(master,LEFT)
        rest=self.ram(master,LEFT)
        uppruta=self.ram(storruta,TOP)
        nerruta=self.ram(storruta,TOP)
        self.atomruta(uppruta,p.atom,p.num)
        if p.down:
            self.stolpe(nerruta)
            self.picture(nerruta,p.down)
            self.ram(nerruta,TOP)
        if p.next:
            self.streck(uppruta)
            self.picture(rest,p.next)

    def show(self,p):
        """ Ritar hela bilden. Parametrar: p (referens till datastrukturen som ska ritas) """
        if self.root!=None:
            self.root.destroy()
        self.root=Tk()
        Label(self.root,text="  ",font=self.stor,bg="white").pack(side=LEFT,fill=Y)
        self.picture(self.root,p)
        #mainloop() #Kommentera bort om du anv. IDLE (IDLE har egen mainloop())
