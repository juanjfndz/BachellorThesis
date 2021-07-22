This is a repository of my final thesis in physics. There are 4 elements

TFG_code.py: 

  It is a work in which a modular code is developed in Python that allows to extract information from the simulations of the EAGLE project. 

In addition, two Jupyter (Colab) Notebooks are shown as examples that allow the work carried out in this study to be visualised.


TFG_Data_Study: 

  This Notebook is a study of the data provided by the simulation databases of the EAGLE project. This gives an idea of what was and what was not to          be developed for the code.



TFG_Notebooks_example: 

  The second Notebook is a complete example of how to use the developed code. It allows to obtain raw data from the EAGLE simulations, work until the necessary information is obtained, save it in a dataset of its own and finally visualise the data obtained. 


To talk in a little more detail about the code: 

This is a code that presents a class and 3 functions. The class "Data_snapnum" and the function "Galaxy_to_past" are the main tools to extract information from the simulations. The other two functions are examples of what can be done once the information has been extracted.

The idea of this code is to be modular. If instead of using the last two functions, you want to develop other functions for other purposes, you don't have to change the whole code. Just create the new functions and rely on the class "Data_snapnum" and "Galaxy_to_past".
