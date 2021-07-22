# -*- coding: utf-8 -*-
"""		TFG_code.py -- version 4.14.0 -- 		"""

import os                                                                       # Librery to work with our system operative.
import h5py                                                                     # Library to work with hdf5 files.
import gc                                                                       # Library to clean trash data.
import numpy as np                                                              # Library to work with data.
import eagleSqlTools as eagle                                                   # Library with EAGLE's tools.

list_sim = ["RefL0012N0188", "RefL0025N0376",
	    "RefL0050N0752", "RecalL0025N0752"]                                # List of EAGLE´s simulations.

sim = list_sim[0]							   # Simulation selected



"""	# Class Data_snapnum	"""



# Commented out IPython magic to ensure Python compatibility.
class Data_snapnum():                                                           # Class to get the data from a snapshot of a simulation.
  def __init__(self, simulation, snapnum, path,				    # Init function.
               user_name="<user_name>", password="<password>"):
    """
      To load the class is requested:
        simulation : The simulation from which the information is to be obtained. 
        snapnum    : The number related to the snapshot of this simulation. 
        Path	  : Path where the .hdf5 file with the snapshot information.
        user_name, password: The username and password of an EAGLE account.
    """

    # The input data is stored in the object.
    self.user    = {'username':user_name, 'password': password}                 # EAGLE User information.
    self.sim     = simulation                                                   # Name of the simulation (e.g. RefL0012N0188).
    self.snapnum = snapnum                                                      # Snapshot number (0, ..., 28).
    self.path    = path                                                         # Path of the database.
    self.subpath = "snap_"+self.path[len("/%s/snapshot"%(self.sim)):]           # Subpath within the path that is used to obtain certain data.
    self.nfiles  = len(os.listdir('%s/'%(path)))				    # nfiles inside of the snapshot.

    with h5py.File('%s/%s.%i.hdf5'%(self.path, self.subpath, 0), 'r') as f:     # Certain interesting constants are loaded.
      self.a     = f['Header'].attrs.get('Time')                                # Scale factor.
      self.h     = f['Header'].attrs.get('HubbleParam')                         # Hubble constan (H_0/100km s^{-1} Mpc^{-1}): 0.667
      self.m_DM  = f['Header'].attrs.get('MassTable')[1]*self.h**[-1]           # DarkMatter mass (M_sol) CGS = 1.989e43 g
    
      # Other Certain interesting constants are loaded (rho_c)
      self.rho_c      = 30000/((3.086e+19)**2*6.6743e-8*8*np.pi)*(self.h)**2    # critical density (g/cm^3)
      Omega_matter    = f['Header'].attrs.get('Omega0')                         # matter density (g/cm^3)
      self.rho_matter = self.rho_c*Omega_matter*self.a**(-3)                    # average matter density at this snapshot.
    
      aexp = f['PartType%i/%s'%(0, 'Coordinates')].attrs.get('aexp-scale-exponent')
      hexp = f['PartType%i/%s'%(0, 'Coordinates')].attrs.get('h-scale-exponent')

      self.boxsize = f['Header'].attrs.get('BoxSize')*self.a**(aexp)*self.h**(hexp)# Size of the universe (Mpc)

      del(aexp, hexp, Omega_matter)                                             # Elimination of time variables
      gc.collect()                                                                

    if self.snapnum != 0:                                                       # load the Galaxies catalogue.
      self.load_catalogue()


  def load_catalogue(self):                                                     # Function to load certain galaxy data.
    """
      Function that allows to load the different galaxies with certain 
      properties such as their GroupNumber and SubGroupNumber. 

      This is valid for all that is different from snapnum 0, because it is 
      from that snapnum that the GroupNumber and SubGroupNumber that the 
      GroupNumber and SubGroupNumber are registered.

      Query information:
        Snapnum: Snapshot's number
        Galaxy ID: Number for each galaxy at snapnum number
        GroupNumber: Número para diferenciar regiones de la simulación
        SubGroupNumber: Número para diferencias subregiones
        KappaCoRot: Disc parameter (if KappaCoRot > 0.4 -> Galaxy Disc-like)
    """

    con   = eagle.connect(user     = self.user['username'], 
                          password = self.user['password'])                      # Connection with the database
    
    query = "SELECT \
                MK.GalaxyID, \
                SH.GroupNumber,\
                SH.SubGroupNumber, \
                SH.Vmax, \
                SH.MassType_Star, \
                SH.Mass, \
		SH.GasSpin_z, \
                MK.KappaCoRot \
          FROM \
                %s_SubHalo AS SH, \
                %s_MorphoKinem AS MK \
          WHERE \
                SH.GalaxyID = MK.GalaxyID AND \
                SH.MassType_Star >= 1E09 AND \
                MK.KappaCoRot >= 0.4 AND \
                SH.SubGroupNumber = 0 AND \
                SH.SnapNum = %i \
          ORDER BY \
                SH.GalaxyID"%(self.sim, self.sim, self.snapnum)                 # Query with the conditions for massive disc galaxies

    self.catalogue = eagle.execute_query(con , query)                           # Catalogue variable

    del(con, query)                                                             # Elimination of time variables
    gc.collect()                                                                

  def read_dataset(self, itype, att):                                           # Function to read the data  
    """
      itype is the type of the particle:
          PartType0 --> Gas particle data
          PartType1 --> Dark matter particle data
          PartType4 --> Star particle data
          PartType5 --> Black hole particle data
      att is the attribute which is looking for.
    """

    with h5py.File('%s/%s.%i.hdf5'%(self.path, self.subpath, 0), 'r') as f:     # File upload to know the number particles.
      n_particles = f['Header'].attrs.get('NumPart_Total')[itype]               # Particles' number.

      if att == 'Coordinates' or att == 'Velocity':
        data = np.ones((n_particles, 3))                                        # Output array for Coordinates and Velocity cases.

      else:
        data = np.ones(n_particles)                                             # Output array rest of cases

      del(n_particles)                                                          # Elimination of time variables
      gc.collect

    if itype==1 and att=='Mass':                                                # DarkMatter Mass case
      data      *= self.m_DM                                                    # Array with all the values of the masses of all these particles
      data.dtype = [('Mass', data.dtype)]                                       # Add the label
      return data

    # The rest
    else:
      try:                                                                      # There are failures when the specified itype is not present, so a it is used a try-except estructure
        count = 0
        for i in range(self.nfiles):                                            # Loop over each file and extract the data.
        
          f      = h5py.File('%s/%s.%i.hdf5'%(self.path, self.subpath, i), 'r') # Load the i-st File
          tmp    = f['PartType%i/%s'%(itype, att)][...]                         # The att in that File
          data[count:count+len(tmp)] = tmp                                      # Saving all the data together
          count += len(tmp)
        
        """
        In the paper (arXiv: 1706.09899, part 4.1) repeat these calcs for each
        nfile (i.e. i value) but they didn't save in any variable, they overwrite all
        these variables, so it's fast if we do only in the last one, or first or whatever.

        I choose the last one.
        """
        aexp = f['PartType%i/%s'%(itype, att)].attrs.get('aexp-scale-exponent')
        hexp = f['PartType%i/%s'%(itype, att)].attrs.get('h-scale-exponent')

        f.close()
        del(f, tmp, count, i)                                                   # Elimination of time variables
        gc.collect()
        
        # convert comovil to movil units and posible h multiplications
        if att != 'ParticleIDs' and data.dtype != np.int32 and data.dtype != np.int64:                 
          data = np.multiply(data, self.a**aexp * self.h**hexp, dtype='f8')

        del(aexp, hexp, tmp)                                                    # Elimination of time variables
        gc.collect()

      except KeyError:                                                          # In case there are no particles of that itype in this snapshot
        if att == 'Coordinates' or att == 'Velocity':                           # Coordinates case
          data = np.ones(3)*np.nan
        else:
          data = np.array([np.nan])
      
      finally:
        if att == 'Coordinates' or att == 'Velocity':                           # Coordinates and velocity case
          data.dtype = [(att+'_%i'%(i), data.dtype) for i in [0, 1, 2]]

        elif att == 'ParticleIDs':                                              # Particle case
          data.dtype = [('ParticleIDs', '<u8')]

        else:                                                                   # Rest
          data.dtype = [(att, data.dtype)]

        return data

  def periodicity(self, array, point, center=False):
    """
    A function that allows particles to be centred around the coordinates 
    of a point or a particle, based on the periodicity property of the universe.
    """

    if point.dtype == np.dtype([('ParticleIDs', '<u8'), ('Coordinates_0', '<f8'), 
                                ('Coordinates_1', '<f8'), ('Coordinates_2', '<f8'), 
                                ('Mass', '<f8'), ('itype', 'i1')]):             # Particles case
      for i in [0, 1, 2]:                                                       # For each space-dimension
        i_point = point['Coordinates_%i'%(i)]
        array['Coordinates_%i'%(i)] -= i_point                                  # Particle as a reference point

        mask = array['Coordinates_%i'%(i)] > self.boxsize/2                     # For all particle beyond L/2...
        array['Coordinates_%i'%(i)] -= mask.astype(np.int)*self.boxsize         # ... is placed on the other side
        del(mask)                                                               # Elimination of time variables
        gc.collect()

        mask = array['Coordinates_%i'%(i)] < -self.boxsize/2                    # Fot all particle beyond -L/2
        mask = mask.astype(np.int)*self.boxsize
        array['Coordinates_%i'%(i)] += mask                                     # ... is placed on the other side
        
        if not(center):                                                         # IF center = False -> the coordinates are retrieved
          array['Coordinates_%i'%(i)] += i_point
        
        del(mask, i_point)                                                               # Elimination of time variables
        gc.collect()

      del(i, point)
      gc.collect()
      return array

    else:                                                                       # Point case (The rest is the same)
      for i in [0, 1, 2]:
        array['Coordinates_%i'%(i)] -= point[i]

        mask = array['Coordinates_%i'%(i)] > self.boxsize/2
        array['Coordinates_%i'%(i)] -= mask.astype(np.int)*self.boxsize
        
        del(mask)
        gc.collect()

        mask = array['Coordinates_%i'%(i)] < -self.boxsize/2
        mask = mask.astype(np.int)*self.boxsize
        array['Coordinates_%i'%(i)] += mask

        if not(center):
          array['Coordinates_%i'%(i)] += point[i]
        
      del(mask, i)
      gc.collect()

      return array

  def particles_prop(self, att=None, itype=None, gn=None, sgn=None):            # Function to obtain the particle data of a certain sector
    """
      A galaxy is defined by its GroupNumber and SubGroupNumber,
      extract the coordinates of all particles of a selected type.
      Coordinates are then wrapped around the centre to account for periodicity.
        * simulation, e.g. RefL0012N0188
        * snapnum, snapshot number (0, ..., 28)
        * itype is the i-st PartType:
            PartType0 --> Gas particle data
            PartType1 --> Dark matter particle data
            PartType4 --> Star particle data
            PartType5 --> Black hole particle data
        * gn is de GroupNumber:
            Tag with a maximum of 2**30 particles and runs 1 to N 
        * sgn is the SubGroupNumber:

        * centre is the centre of the preiodicity
    """
    if itype == None and att == None:
      with h5py.File('%s/%s.%i.hdf5'%(snap_0.path, snap_0.subpath, 0), 'r') as f:# File upload to conover the number particles.
        n_particles = f['Header'].attrs.get('NumPart_Total')                    # Particles' number
        DF          = np.ones(np.sum(n_particles), 
	      		      dtype=np.dtype([('ParticleIDs', 'u8'), ('Coordinates_0', '<f8'), \
                                              ('Coordinates_1', '<f8'), ('Coordinates_2', '<f8'), \
                                              ('Mass', '<f8'), ('itype', 'i1')]))# Output array vacuum 
	
      count = 0
      for itype, i in zip([0, 1, 4, 5], [0, 1, 2, 3]):		             # For each itype

        if n_particles[i] == 0:
          continue

        for att in ['ParticleIDs', 'Coordinates', 'Mass', 'itype']:		     # For each att add in DF the data in the right position
          if att == 'Coordinates':
            data = snap_0.read_dataset(itype, att)
            DF['Coordinates_0'][count:count+n_particles[i]] = data['Coordinates_0'][..., 0]
            DF['Coordinates_1'][count:count+n_particles[i]] = data['Coordinates_1'][..., 0]
            DF['Coordinates_2'][count:count+n_particles[i]] = data['Coordinates_2'][..., 0]
            del(data)
            gc.collect()

          elif att == 'itype':
            DF['itype'][count:count+n_particles[i]] *= itype

          else:
            DF[att][count:count+n_particles[i]] = snap_0.read_dataset(itype= itype,att= att)

        count += n_particles[i]

      del(count, n_particles, att, itype, i)                                      # Elimination time variables
      gc.collect()


    elif itype == None:
      dicc_dtype = {'ParticleIDs':[('ParticleIDs'  , 'u8')] , 
                    'Coordinates':[('Coordinates_0', '<f8') ,
                                   ('Coordinates_1', '<f8') , 
                                   ('Coordinates_2', '<f8')],
                    'Mass'       :[('Mass'         , '<f8')],
                    'itype'      :[('itype'        , 'i1')]}


      with h5py.File('%s/%s.%i.hdf5'%(snap_0.path, snap_0.subpath, 0), 'r') as f: # File upload to conover the number particles.
        n_particles = f['Header'].attrs.get('NumPart_Total')                      # Particles' number
        DF          = np.ones(np.sum(n_particles), dtype=np.dtype(dicc_dtype[att]))# Output array vacuum 
	
      count = 0
      for itype, i in zip([0, 1, 4, 5], [0, 1, 2, 3]):			      # For each itype

        if n_particles[i] == 0:
          continue

        if att == 'Coordinates':
          data = snap_0.read_dataset(itype, att)
          DF['Coordinates_0'][count:count+n_particles[i]] = data['Coordinates_0'][..., 0]
          DF['Coordinates_1'][count:count+n_particles[i]] = data['Coordinates_1'][..., 0]
          DF['Coordinates_2'][count:count+n_particles[i]] = data['Coordinates_2'][..., 0]
          del(data)
          gc.collect()

        elif att == 'itype':
          DF['itype'][count:count+n_particles[i]] *= itype

        else:
          DF[att][count:count+n_particles[i]] = snap_0.read_dataset(itype= itype,att= att)

        count += n_particles[i]

      del(count, n_particles, itype, i)                                           # Elimination time variables
      gc.collect()


    else:
      with h5py.File('%s/%s.%i.hdf5'%(snap_0.path, snap_0.subpath, 0), 'r') as f: # File upload to conover the number particles.
        n_particles = f['Header'].attrs.get('NumPart_Total')[itype]               # Particles' number
        DF          = np.ones(n_particles, dtype=np.dtype([('ParticleIDs', 'u8'), ('Coordinates_0', '<f8'), \
                                                           ('Coordinates_1', '<f8'), ('Coordinates_2', '<f8'), \
                                                           ('Mass', '<f8'), ('itype', 'i1')]))    # Output array
	
      count = 0
      if n_particles[i] == 0:
        pass
      else:
        for att in ['ParticleIDs', 'Coordinates', 'Mass', 'itype']:
          if att == 'Coordinates':
            data = snap_0.read_dataset(itype, att)
            DF['Coordinates_0'][count:count+n_particles[i]] = data['Coordinates_0'][..., 0]
            DF['Coordinates_1'][count:count+n_particles[i]] = data['Coordinates_1'][..., 0]
            DF['Coordinates_2'][count:count+n_particles[i]] = data['Coordinates_2'][..., 0]
            del(data)
            gc.collect()

          elif att == 'itype':
            DF['itype'][count:count+n_particles[i]] *= itype

          else:
            DF[att][count:count+n_particles[i]] = snap_0.read_dataset(itype= itype,att= att)

          count += n_particles[i]

        del(count, n_particles, att)                                            # Elimination time variables
        gc.collect()


    if gn != None and sgn != None and self.snapnum != 0:                        # There is the possibility to add a mask over the GroupNumber and SubGroupNumber.

      gns  = self.read_dataset(itype, 'GroupNumber')['GroupNumber']
      sgns = self.read_dataset(itype, 'SubGroupNumber')['SubGroupNumber']
      mask = np.logical_and(gns == gn, sgns == sgn)                             # Mask
      del(gns, sgns)
      gc.collect()

      DF = DF[mask]
      del(mask)
      gc.collect()

      DF.dtype=np.dtype([('ParticleIDs', 'u8'), ('Coordinates_0', '<f8'), \
                         ('Coordinates_1', '<f8'), ('Coordinates_2', '<f8'), \
                         ('Mass', '<f8'), ('itype', 'i1')])          

    return DF



"""	# Galaxy_to_past		"""



def Galaxy_to_past(GalaxyID, snap_1, snap_2):                                   # Function to obtain center and the radio of the sphere with all the particles at snap_2 of the galaxy at snap_1
  """
    Function that allows to study the overdensity formed by the particles at 
    snapnum 0 that at snapnum 28 form a galaxy.

    snap_1 --> Galaxy snapshot
    snap_2 --> Past snapshot
    GalaxyID --> ID of the galaxy at snap_1
    
    We only consider to use the darkmatter 
  """


  mask_gn_sgn = snap_1.catalogue['GalaxyID'] == GalaxyID                        # Galaxy Mask
  gn, sgn     = snap_1.catalogue[['GroupNumber', 
                                  'SubGroupNumber']][mask_gn_sgn][0]            # GroupNumber and SubGroupNumber of the Galaxy

  del(mask_gn_sgn)                                                              # Elimination time variables
  gc.collect()
  
  gns    = snap_1.read_dataset(itype= 1, att='GroupNumber')['GroupNumber']
  sgns   = snap_1.read_dataset(itype= 1, att='SubGroupNumber')['SubGroupNumber']
  mask_1 = np.logical_and(gns == gn, sgns == sgn)

  n_particles   = np.sum(mask_1)                    
  ParticleIDs_1 = snap_1.read_dataset(itype= 1, att='ParticleIDs')['ParticleIDs'][mask_1] # Array with the IDs of the galaxy at snap_1
    
  del(gns, sgns, mask_1)                                                        # Elimination time variables
  gc.collect()

  ParticleIDs_2 = snap_2.read_dataset(itype = 1, att='ParticleIDs')['ParticleIDs'] # All IDs Particles at snap_2
  mask_2        = np.in1d(ar1 = ParticleIDs_2, ar2 = ParticleIDs_1)             # Mask to looking where are the galaxy's particles on snapshot 2.

  del(ParticleIDs_1, ParticleIDs_2)                                             # Elimination time variables
  gc.collect()
  
  Coordinates = snap_2.read_dataset(itype = 1, att='Coordinates')[mask_2][:,0]  # Coordinates of the galaxys particles at snap_2
  del(mask_2)
  gc.collect()
  Coordinates = snap_2.periodicity(Coordinates, Coordinates[0].copy())          # Periodicity around first particle.    

  Mass_T      = snap_2.m_DM*n_particles				            # Its mass
  del(n_particles)
  gc.collect()

  Mass_center = np.array([np.sum(Coordinates['Coordinates_%i'%(i)]*snap_2.m_DM) for i in [0, 1, 2]])/Mass_T# Mass_center
  del(Mass_T)
  gc.collect()
  
  Radios = np.sum([(Coordinates['Coordinates_%i'%(i)]-Mass_center[i])**2 for i in [0, 1, 2]], axis=0)  # the radios of the galaxy's particle with the new mass center
  del(Coordinates)
  gc.collect()
  Radio = max(np.sort(Radios)[:int(len(Radios)*0.9)]) 
  del(Radios)
  gc.collect()
  return Radio, Mass_center
  


"""	# Overrho	"""


def Overrho(snap, Radio, center):
  """
	Function to obtain the over-density of the sphere given as input 
        (radio and center) at the snapshot also selected (snap)
  """

  with h5py.File('%s/%s.%i.hdf5'%(snap_0.path, snap_0.subpath, 0), 'r') as f:   # File upload to know the number particles.
    n_particles = f['Header'].attrs.get('NumPart_Total')                        # Particles' number
    Radios_2    = np.ones(np.sum(n_particles), dtype=np.dtype('<f8'))           # Output array vacuum, 

  count = 0
  for i in [0, 1]:                                                              # radius of gas and DM particles with respect the center of mass.
    Coordinates = snap.read_dataset(itype= i, att='Coordinates')
    Coordinates = snap.periodicity(Coordinates, center)
    Radios_2[count:count+n_particles[i]] = np.sum([(Coordinates['Coordinates_%i'%(i)]-center[i])**2 for i in [0, 1, 2]], axis=0)[:,0]

    del(Coordinates)
    gc.collect()
      
    count += n_particles[i]
  
  del(i)
  gc.collect()

  Mass = np.sum(snap.particles_prop(att='Mass')['Mass'][Radios_2 <= Radio])     # Total mass inside of the sphere
  rho  = Mass*1.989e43/((4/3)*np.pi*((Radio)**(1/2)*3.08e24)**3)                # Density of that sphere
  overrho = (rho - snap.rho_matter)/snap.rho_matter                             # Overdensity  
  return overrho, rho, Mass



"""	# AngularMoment      """


def AngularMoment(snap, Radio, center):
  """
	Function to obtain the total vector angular momentum of the sphere
	given as input (radio and center) at the snapshot also selected (snap)
  """
  with h5py.File('%s/%s.%i.hdf5'%(snap_0.path, snap_0.subpath, 0), 'r') as f:   # File upload to know the number particles.
    n_particles = f['Header'].attrs.get('NumPart_Total')                        # Particles' number
    Radios_2    = np.ones(np.sum(n_particles), dtype=np.dtype('<f8'))           # Output array vacuum 

  count = 0
  for i in [0, 1]:                                                              # gas and DM particles radios from the center of mass.
    Coordinates = snap.read_dataset(itype= i, att='Coordinates')		    # radius of gas and DM particles with respect the center of mass.
    Coordinates = snap.periodicity(Coordinates, center)
    Radios_2[count:count+n_particles[i]] = np.sum([(Coordinates['Coordinates_%i'%(i)]-center[i])**2 for i in [0, 1, 2]], axis=0)[:,0]

    del(Coordinates)
    gc.collect()
      
    count += n_particles[i]
  
  del(i)
  gc.collect()
  
  count = 0
  count_angular = 0
  angulars = np.ones((np.sum(Radios_2 <= Radio), 3), dtype=np.dtype('<f8'))
  for i in [0, 1]:                                                              # Calculus the vectorial product between the coordinates and the velocity of the particles inside the region
    mask = [Radios_2[count:count+n_particles[i]] <= Radio][0]
    Coordinates = snap.read_dataset(itype= i, att='Coordinates')[..., 0][mask]
    Coordinates = snap.periodicity(Coordinates, center)
    Velocity    = snap.read_dataset(itype= i, att='Velocity')[..., 0][mask]

    angulars[count_angular:count_angular + np.sum(mask)][..., 0] = Coordinates['Coordinates_1']*Velocity['Velocity_2'] - Coordinates['Coordinates_2']*Velocity['Velocity_1']
    angulars[count_angular:count_angular + np.sum(mask)][..., 1] = Coordinates['Coordinates_2']*Velocity['Velocity_0'] - Coordinates['Coordinates_0']*Velocity['Velocity_2']
    angulars[count_angular:count_angular + np.sum(mask)][..., 2] = Coordinates['Coordinates_0']*Velocity['Velocity_1'] - Coordinates['Coordinates_1']*Velocity['Velocity_0']
   
    del(Coordinates, Velocity)
    gc.collect()
      
    count += n_particles[i]
    count_angular += np.sum(mask)

    if i == 1:								     # times DM mass for DM particles
      angulars[count_angular:count_angular + np.sum(mask)] *= snap.m_DM
    elif i == 0:								     # times gas mass for gas particles
      mass = snap.read_dataset(itype= i, att='Mass')['Mass'][mask][0]
      angulars[count_angular:count_angular + np.sum(mask)] *= mass
      del(mass)
    
    del(mask)

  del(i, Radios_2)
  gc.collect()
  
  angular = np.array([np.sum(angulars[i]) for i in [0, 1, 2]])
  return angular




	
"""# Example storing the data"""

snap_0  = Data_snapnum(simulation=sim, snapnum=0,                        
                      path='%s/snapshot_000_z020p000'%(sim))           	        # first snapshot load
snap_28 = Data_snapnum(simulation=sim, snapnum=28, 
                      path='%s/snapshot_028_z000p000'%(sim))                        # second snapshot load

n_galaxies = len(snap_28.catalogue)
print("%s: %s disc massive galaxies: "%(sim, n_galaxies))

h5f = h5py.File('TFG_%s.h5'%(sim), 'a')				                # Creation of the h5 file for all the information needed
try:									        # Creation of the simulation section
  IDf = h5f.create_dataset(sim, data=snap_28.catalogue['GalaxyID'])
except:
  IDf = h5f['%s'%(sim)]


IDf.attrs['Mass']       = snap_28.catalogue['Mass']
IDf.attrs['MassStar']   = snap_28.catalogue['MassType_Star']
IDf.attrs['Vmax']       = snap_28.catalogue['Vmax']
IDf.attrs['GasSpin_z']  = snap_28.catalogue['GasSpin_z']
IDf.attrs['KappaCoRot'] = snap_28.catalogue['KappaCoRot']

                      
print("Save 1st Info  --> Done")

radio       = np.zeros(n_galaxies, dtype= np.dtype('<f8'))
mass_center = np.zeros((n_galaxies, 3), dtype= np.dtype('<f8'))
Mass        = np.zeros(n_galaxies, dtype= np.dtype('<f8'))
rho         = np.zeros(n_galaxies, dtype= np.dtype('<f8'))
overrho     = np.zeros(n_galaxies, dtype= np.dtype('<f8'))

for i in range(n_galaxies):						        # Radio, Mass_center and overrho
  radio[i], mass_center[i]    = Galaxy_to_past(GalaxyID = snap_28.catalogue['GalaxyID'][i], 
		                               snap_1   = snap_28, snap_2 = snap_0) 
  overrho[i], rho[i], Mass[i] = Overrho(snap = snap_0, Radio = radio[i], center = mass_center[i])# overrho

IDf.attrs['radio_sphere']   = radio
IDf.attrs['center_mass']    = mass_center
IDf.attrs['Snap_0_mass']    = Mass
IDf.attrs['Snap_0_rho']     = rho
IDf.attrs['Snap_0_overrho'] = overrho

print("Save 2nd Info  --> Done\n")

atrib = lambda att: IDf.attrs.get(att)
angular = np.zeros((n_galaxies,3), dtype= np.dtype('<f8'))

for i in range(n_galaxies):
  angular[i] = AngularMoment(snap = snap_0, Radio = atrib('radio_sphere')[i], center = atrib('center_mass')[i])

IDf.attrs['Snap_0_angular'] = angular
print("\nSave 3rd Info  --> Done")

h5f.close()

print(input('\n Finised.s'))






