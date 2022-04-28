#!/usr/bin/env python3
import sys
import matplotlib
matplotlib.use('WebAgg')
from matplotlib import pyplot as plt

seqs='''
>name1
fffffgghhg
>name2
fffffggggg
>name3
alfdjsvfdv
'''

separacion=0.3

def read_fasta(txt):
    '''
    Reads sequences in fasta format and returns them in a dictionary.
    Takes a string of sequences in FASTA format and returns a dictionary of the form {"name1":"sequence1","name2",sequence2"...}
    Parameters:
    txt (str): Sequences

    Returns:
    seq_dict (dict): a dictionary of the form {"name1":"sequence1","name2",sequence2"...}                                                       
    '''
    seq_dict={}
    j=2
    es_nombre=True
    array=[]
    nombre=''
    secuencia=''
    while j<len(txt):
        letter=txt[j]
        if not es_nombre:
            if letter=='\n':
                secuencia=secuencia
            elif letter=='>':
                #print('empieza a ser nombre')
                seq_dict[nombre]=secuencia
                es_nombre=True  
                nombre=''

            else:
                #print('sigue siendo seq')
                secuencia=secuencia+letter
        else:
            if letter=='\n':
                #print('empieza a ser seq')
                es_nombre=False
                array.append([secuencia,nombre])
                secuencia=''
            else:
                #print('sigue siendo nombre')
                nombre=nombre+letter
        #print(j)
        j=j+1
    seq_dict[nombre]=secuencia
    return seq_dict
#print(seqs)
#A=read_fasta(seqs)
#print(A)

class especie:
  '''a general class of species, each instance of this class will have a position in the tree'''
  def __init__(self):
    pass

  @staticmethod
  def id_actuales(especie1,especie2):
    '''

      Parameters
      ----------
      especie1 : (especie_actual).
      especie2 : (especie_actual). 

      Returns
      -------
      (float)
          %matches between the sequences.

      '''
    secuencia1=especie1.secuencia
    secuencia2=especie2.secuencia
    if len(secuencia1)!=len(secuencia2):
      print('las secuencias no tienen la misma longitud')
      return None
    else:
      n_id=0
      for i in range(len(secuencia1)):
        if secuencia1[i]==secuencia2[i]:
          if secuencia1[i]!='-':
            n_id=n_id+1
      return 100*n_id/len(secuencia1)
  
  @classmethod
  def id_actual1(cls,especie1,especie2):
    '''
      

      Parameters
      ----------
      cls : (class).
      especie1 : (especie_actual).
      especie2 : (ancestro).

      Returns
      -------
      id : (float)
          Estimated identity percentage between especie1 and especie2.

      '''
    secuencia1=especie1.secuencia
    hijo1=especie2.hijo1
    hijo2=especie2.hijo2
    #id_1
    if type(hijo1)==ancestro:
      id_1=cls.id_actual1(especie1,hijo1)
    elif type(hijo1)==especie_actual:
      id_1=cls.id_actuales(especie1,hijo1)
    #id_2
    if type(hijo2)==ancestro:
      id_2=cls.id_actual1(especie1,hijo2)
    elif type(hijo2)==especie_actual:
      id_2=cls.id_actuales(especie1,hijo2)
    id=(id_1+id_2)/2
    return id

  @classmethod
  def id_actual2(cls,especie1,especie2):
    '''
      

      Parameters
      ----------
      cls : (class).
      especie1 : (ancestro).
      especie2 : (especie_actual).

      Returns
      -------
      id : (float)
          Estimated identity percentage between especie1 and especie2.

      '''
    id=cls.id_actual1(especie2,especie1)
    return id
  @classmethod
  def id_ancestros(cls,especie1,especie2):
    '''
      

      Parameters
      ----------
      cls : (class).
      especie1 : (ancestro).
      especie2 : (ancestro).

      Returns
      -------
      id : (float)
          Estimated identity percentage between especie1 and especie2.

      '''
    hijo11=especie1.hijo1
    hijo12=especie1.hijo2
    hijo21=especie2.hijo1
    hijo22=especie2.hijo2
    #id_11 hijo 1 especie 1 vs hijo 1 especie 2
    if type(hijo11)==ancestro and type(hijo21)==ancestro:
      id_11=cls.id_ancestros(hijo11,hijo21)
    elif type(hijo11)==especie_actual and type(hijo21)==ancestro:
      id_11=cls.id_actual1(hijo11,hijo21)
    elif type(hijo11)==ancestro and type(hijo21)==especie_actual:
      id_11=cls.id_actual2(hijo11,hijo21)
    elif type(hijo11)==especie_actual and type(hijo21)==especie_actual:
      id_11=cls.id_actuales(hijo11,hijo21)
    #id_12 hijo 1 especie 1 vs hijo 2 especie 2  
    if type(hijo11)==ancestro and type(hijo22)==ancestro:
      id_12=cls.id_ancestros(hijo11,hijo22)
    elif type(hijo11)==especie_actual and type(hijo22)==ancestro:
      id_12=cls.id_actual1(hijo11,hijo22)
    elif type(hijo11)==ancestro and type(hijo22)==especie_actual:
      id_12=cls.id_actual2(hijo11,hijo22)
    elif type(hijo11)==especie_actual and type(hijo22)==especie_actual:
      id_12=cls.id_actuales(hijo11,hijo22)
    #id_22 hijo 2 especie 1 vs hijo 2 especie 2  
    if type(hijo12)==ancestro and type(hijo22)==ancestro:
      id_22=cls.id_ancestros(hijo12,hijo22)
    elif type(hijo12)==especie_actual and type(hijo22)==ancestro:
      id_22=cls.id_actual1(hijo12,hijo22)
    elif type(hijo12)==ancestro and type(hijo22)==especie_actual:
      id_22=cls.id_actual2(hijo12,hijo22)
    elif type(hijo12)==especie_actual and type(hijo22)==especie_actual:
      id_22=cls.id_actuales(hijo12,hijo22)
    #id_21 hijo 2 especie 1 vs hijo 1 especie 2  
    if type(hijo12)==ancestro and type(hijo21)==ancestro:
      id_21=cls.id_ancestros(hijo12,hijo21)
    elif type(hijo12)==especie_actual and type(hijo21)==ancestro:
      id_21=cls.id_actual1(hijo12,hijo21)
    elif type(hijo12)==ancestro and type(hijo21)==especie_actual:
      id_21=cls.id_actual2(hijo12,hijo21)
    elif type(hijo12)==especie_actual and type(hijo21)==especie_actual:
      id_21=cls.id_actuales(hijo12,hijo21)
    return (id_11+id_12+id_22+id_21)/4
    
  @classmethod
  def identidad(cls,especie1,especie2):
    '''
       

       Parameters
       ----------
       cls : (class).
       especie1 : (especie).
       especie2 : (especie).

       Returns
       -------
       id : (float)
           Estimated identity percentage between especie1 and especie2.
    '''
    if type(especie1)==ancestro and type(especie2)==ancestro:
      id=cls.id_ancestros(especie1,especie2)
    elif type(especie1)==especie_actual and type(especie2)==ancestro:
      id=cls.id_actual1(especie1,especie2)
    elif type(especie1)==ancestro and type(especie2)==especie_actual:
      id=cls.id_actual2(especie1,especie2)
    elif type(especie1)==especie_actual and type(especie2)==especie_actual:
      id=cls.id_actuales(especie1,especie2)
    return id

  @classmethod
  def tevo(cls,especie1,especie2):
      '''

      Parameters
      ----------
      cls : (class).
      especie1 : (especie).
      especie2 : (especie).
      
      Returns
      -------
      (float)
          Value associated with estimated evolutionary time between especie1 and especie2 to graph the x axis of the tree.

      '''
      return 1/cls.identidad(especie1,especie2)
  
  @classmethod
  def matriz_identidad(cls,array_especies):
    '''
      

      Parameters
      ----------
      cls : (class).
      array_especies : (list)
          A list of especie type objects.

      Returns
      -------
      (list)
          A matrix of %identitys whose components M[i][j] are the %identity between the species array_especies[i] and array_especies[j].

      '''
    vector_id=[]
    vector_esp=[]
    i=0
    while i<len(array_especies):
      lista_names=[]
      lista_print=[]
      j=0
      while j<i:
        vector_esp.append([array_especies[i],array_especies[j]])
        vector_id.append(cls.identidad(array_especies[i],array_especies[j]))
        lista_names.append([array_especies[i].name,array_especies[j].name])
        lista_print.append(int(cls.identidad(array_especies[i],array_especies[j])))
        j=j+1
      while j<len(array_especies):
        lista_names.append('')
        lista_print.append('')
        j=j+1
      #print(lista_names)
      #print(lista_print)
      i=i+1
    #print(vector_id)
    compuesto_ordenado=sorted(zip(vector_id,vector_esp), key=lambda pair: pair[0])
    vector_id_ordenado=[id for (id,esp) in compuesto_ordenado]
    vector_esp_ordenado=[par_esp for id,par_esp in compuesto_ordenado]
    return vector_esp_ordenado[::-1]

  

  @classmethod
  def compress_step(cls,array_especies,graficar=False):
    '''
      

      Parameters
      ----------
      cls : (class).
      array_especies : (list)
          A list of especie type objects.
      graficar : Boolean, optional
          Is set to True the second time this command is run, when the modern species are already ordered acording to phylogeny, and assigns positions in the graph to all ancesters. The default is False.

      Returns
      -------
      nuevo_array : (list)
          New array of especie type objects which contains the ancester of species which were closer to each other than to any other in the previous array, those species who dont have a match in this step are added to the array alone.

      '''
    vector_ordenado=cls.matriz_identidad(array_especies)
    nuevo_array=[]
    especies_usadas=[]
    for par in vector_ordenado:
      if par[0] not in especies_usadas and par[1] not in especies_usadas:
        nuevo_array.append(ancestro(par[0],par[1]))
        especies_usadas.append(par[0])
        especies_usadas.append(par[1])
        if graficar==True:
          if type(par[0])==ancestro and type(par[1])==ancestro:
            par[0].asignar_posicionY()
            par[1].asignar_posicionY()
          elif type(par[0])==especie_actual and type(par[1])==ancestro:
            par[1].asignar_posicionY()
          elif type(par[1])==especie_actual and type(par[0])==ancestro:
            par[0].asignar_posicionY()
          A=ancestro(par[0],par[1])
          A.asignar_posicionY()
          ax.plot([par[0].posicionX,A.posicionX],[par[0].posicionY,A.posicionY],'b')
          ax.plot([par[1].posicionX,A.posicionX],[par[1].posicionY,A.posicionY],'b')
      elif par[0] not in especies_usadas:
        nuevo_array.append(par[0])
        especies_usadas.append(par[0])
      elif par[1] not in especies_usadas:
        nuevo_array.append(par[1])
        especies_usadas.append(par[1])
      usadas_name=[]
      for esp in especies_usadas:
        usadas_name.append(esp.name)
    #print('array=',[x.name for x in nuevo_array])
    return nuevo_array

  @classmethod
  def root_tree(cls,array):
    '''
      

      Parameters
      ----------
      cls : (class).
      array :(list)
          The starting array, consisting of only especie_actual type objects.

      Returns
      -------
      array : (list)
          A list containing only one object, an ancestro type object which is the ancester of all the starting species.

      '''
    while len(array)>1:
      array=cls.compress_step(array)
    return array
    
  @staticmethod
  def tiene_ancestros(array):
      '''
      

      Parameters
      ----------
      array : (list)
          An array of especie type objects.

      Returns
      -------
      bool
          If the array contains ancestro type objects returns True.

      '''
      for elemento in array:
        if type(elemento)==ancestro:
          return True
      return False

  @classmethod
  def flat(cls,array):
    '''
      

      Parameters
      ----------
      cls : (class).
      array : (list)
          An array of especie type objects.

      Returns
      -------
      array : (list)
          An array of especie type objects, with the two "children" of all the ancestro type objects that were in the imput array, recursibly flats the array until no more ancestro type objects remain in the array, at which point the array is ordered.he ancestro type objects that were in the imput array, recursibly flats the array until no more ancestro type objects remain in the array, at which point the array is ordered.

      '''
    while cls.tiene_ancestros(array):
      array_nuevo=[]
      for elemento in array:
        if type(elemento)==ancestro:
          array_nuevo.append(elemento.hijo1)
          array_nuevo.append(elemento.hijo2)
        else:
          array_nuevo.append(elemento)
      array=array_nuevo
    return array

  @classmethod
  def asignar_posiciones(cls,array_ordenado):
    '''
      Asigns positions in a cartesian axis to each especie_actual instances.

      Parameters
      ----------
      cls : (class).
      array_ordenado : (list)
          An ordered array of especie_actual type objects in which any close especies in indexes are close evolutionarily so that the tree doesnt intersect itself.

      Returns
      -------
      None.

      '''
    for i in range(len(array_ordenado)):
      array_ordenado[i].posicionY=i*separacion

  @classmethod
  def graf_tree(cls,array_ordenado):
    '''
      Runs Compress step until it reaches the root but assigning coordinates to each especie instance.

      Parameters
      ----------
      cls : (class).
      array_ordenado :(list)
          An ordered array of especie_actual type objects in which any close especies in indexes are close evolutionarily so that the tree doesnt intersect itself.

      Returns
      -------
      None.

      '''
    cls.asignar_posiciones(array_ordenado)
    
    #print('entrada')

    while len(array_ordenado)>1:
      array_ordenado=cls.compress_step(array_ordenado,True)  


class especie_actual(especie):
  ''' An especie type object with a known sequence.'''
  
  def __init__(self,secuencia,name):
    '''
    Asigns the values given as an imput at the moment of creating the instance to the instance itself.

    Parameters
    ----------
    secuencia : (str).
    name : (str).

     Returns
     -------
    None.

    '''
    self.name=name
    self.secuencia=secuencia
    self.posicionX=0
class ancestro(especie):
  ''' An especie type object which is the ancester of two other especie type objects.
  '''
  def __init__(self,hijo1,hijo2):
    '''
      Asign the values given as an imput at the moment of creating the instance to the instance itself.

      Parameters
      ----------
      hijo1 : (especie).
      hijo2 : (especie).

      Returns
      -------
      None.

      '''
    self.hijo1=hijo1
    self.hijo2=hijo2
    self.edad=self.tevo(hijo1,hijo2)
    self.name='A('+hijo1.name+','+hijo2.name+')'
    self.posicionX=-self.edad
  def asignar_posicionY(self):
    '''
    Asigns the position of the ancester depending on the position of its """children""".

    Returns
    -------
    None.

    '''
    self.posicionY=(self.hijo1.posicionY+self.hijo2.posicionY)/2

fig,ax=plt.subplots()
def make_tree(seq_dict):
    '''
    Plot a phylogenetic tree in a matplotlib.

    Parameters
    ----------
    seq_dict : dict
       a dictionary of the form {"""name1""":"""sequence1""","""name2""","""sequence"""...}.

    Returns
    -------
    None.

    '''
    array_especies=[]
    for name,seq in dict.items(seq_dict):
            array_especies.append(especie_actual(seq,name))
    raiz=especie.root_tree(array_especies)
    array_ordenado=especie.flat(raiz)
    especie.graf_tree(array_ordenado)
    for elemento in array_ordenado:
        #print(elemento.name+'=',elemento.secuencia)
        ax.text(elemento.posicionX,elemento.posicionY,elemento.name)
#make_tree(A)

#b='''
#>seq1
#AAC
#>seq2
#ATC
#>seq3
#AGG
#'''
#B=read_fasta(b)
#make_tree(B)

with open(sys.argv[1]) as f:
    texto=f.read()
    seq_dict=read_fasta(texto)
    make_tree(seq_dict)
plt.savefig(sys.argv[2])

