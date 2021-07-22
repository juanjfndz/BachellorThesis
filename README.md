This is a repository of my final thesis in physics. There are 4 elements

TFG_code.py: 

  This is the main code. It is a library that currently has: 1 class and 3 functions.

  These tools are the main part of my final degree work.  The aim of these is to extract and manipulate information from a dataset. The dataset information is 
  provided by the EAGLE project. This is a state-of-the-art project on simulations of galaxy evolution and formation.

  The code can be classified into two parts: 

    The "Data_snapnum" class and the "Galaxy_to_past" function: These are the main tools of this code. The class allows to connect to the EAGLE dataset and extract
    information. The possible outputs of "Data_snapnum" are focused on the possible needs of the rest of the functions. "Galaxy_to_past" function allows to connect
    information at different times of the simulation. Knowing the position of certain particles at a certain time.

    The "Overrho" and "AngularMoment" functions: These are two examples of what can be done with this code. With the previous tools developed, it only remains to 
    see what you want to investigate. This is the advantage of this code, it is modular. Being modular allows you to create new functions without having to modify 
    the code, just add more functions.



EAGLERawData_Study.ipynb: 

  This is the notebook before the creation of TFG_code.py. This is a notebook focused on the study of raw data. With this notebook you can study how are the outputs  
  of the dataset provided by EAGLE. In addition, it allows us to see what is not in the EAGLE dataset and to set as a goal to obtain it.
  
  

TFG_code_example: 

  This second notebook is an example of how you would work with the TFG_code.py library.  This shows that a very elaborate study can be carried out in just a few 
  lines. Moreover, as it is modular, new tools could be created based on those already developed.

  In addition to using the library, the results obtained in the final degree project are also shown.



TFG_alu0101067766.pdf:

  Finally, the working memory has also been uploaded. In this memory, the meaning of what is asked for in the code, how it can be obtained and what the results   
  obtained imply. 
