# -*- coding: utf-8 -*-
"""		TFG_code.py -- version 4.14.0 -- 		"""

# Libreries
import os                                                                       
import h5py                                                                     
import gc                                                                       
import numpy as np                                                              
import eagleSqlTools as eagle    

"""	# Class Data_snapnum	"""

# Class to get the data from a snapshot of a simulation.
class Data_snapnum():                                                           
  def __init__(self, simulation, snapnum, path,				   
               user_name="<user_name>", password="<password>"):
    """
      To load the class is requested:
        simulation : The EAGLE´s simulation
        snapnum    : The number related to the snapshot of this simulation. 
        Path	   : Path where the .hdf5 file with the snapshot information.
        user_name, password: The username and password of an EAGLE account.
    """

    # The input data is stored in the object.
    self.user    = {'username':user_name, 'password': password}                 
    self.sim     = simulation                                                   
    self.snapnum = snapnum                                                      
    self.path    = path                                                         
    self.subpath = "snap_"+self.path[len("/%s/snapshot"%(self.sim)):]           
    self.nfiles  = len(os.listdir('%s/'%(path)))				    

    # Certain interesting constants are loaded.
    with h5py.File('%s/%s.%i.hdf5'%(self.path, self.subpath, 0), 'r') as f:  
      """
      # Scale factor
      # Hubble constan: 0.667
      # DarkMatter mass (M_*)
      """

      self.a     = f['Header'].attrs.get('Time')          
      self.h     = f['Header'].attrs.get('HubbleParam')               
      self.m_DM  = f['Header'].attrs.get('MassTable')[1]*self.h**[-1] 
    
      # Other Certain interesting constants are loaded (rho_c)
      """
      # critical density (g/cm^3)
      # matter density   (g/cm^3)
      # average matter density at this snapshot.
      """

      self.rho_c      = 30000/((3.086e+19)**2*6.6743e-8*8*np.pi)*(self.h)**2    
      Omega_matter    = f['Header'].attrs.get('Omega0')                         
      self.rho_matter = self.rho_c*Omega_matter*self.a**(-3)                    
    
      aexp = f['PartType%i/%s'%(0, 'Coordinates')].attrs.get('aexp-scale-exponent')
      hexp = f['PartType%i/%s'%(0, 'Coordinates')].attrs.get('h-scale-exponent')

      # Size of the universe (Mpc)
      self.boxsize = f['Header'].attrs.get('BoxSize')*self.a**(aexp)*self.h**(hexp)

      del(aexp, hexp, Omega_matter)                                             
      gc.collect()                                                                

    # load the Galaxies catalogue. (Snapshot 0 hasn´t)
    if self.snapnum != 0:                                                       
      self.load_catalogue()

  # Function to load certain galaxy data.
  def load_catalogue(self):                                                     
    """
      Function that allows to load the different galaxies with certain 
      properties such as their GroupNumber and SubGroupNumber. 

      This is valid for all that is different from snapnum 0, because it is 
      from that snapnum that the GroupNumber and SubGroupNumber that the 
      GroupNumber and SubGroupNumber are registered.

      Query information:
        Snapnum:         Snapshot's number
        Galaxy ID:       Number for each galaxy at snapnum number
        GroupNumber:     Halo's number
        SubGroupNumber:  Galaxy's number
        KappaCoRot: Disc parameter (if KappaCoRot > 0.4 -> Galaxy Disc-like)
    """

    # Connection with the database
    con   = eagle.connect(user     = self.user['username'], 
                          password = self.user['password'])                      
    # Query with the conditions for massive disc galaxies
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
                SH.GalaxyID"%(self.sim, self.sim, self.snapnum)                 
    # Our Catalogue
    self.catalogue = eagle.execute_query(con , query)                           

    del(con, query)                                                             
    gc.collect()                                                                

  # Function to read the data  
  def read_dataset(self, itype, att):                                          
    """
      itype is the type of the particle:
          PartType0 --> Gas particle data
          PartType1 --> Dark matter particle data
          PartType4 --> Star particle data
          PartType5 --> Black hole particle data
      att is the attribute which is looking for.
    """
    
    # File upload to know the number particles.
    with h5py.File('%s/%s.%i.hdf5'%(self.path, self.subpath, 0), 'r') as f:  
      # Particles' number.   
      n_particles = f['Header'].attrs.get('NumPart_Total')[itype]               

      # Output array for Coordinates and Velocity cases.
      if att == 'Coordinates' or att == 'Velocity':
        data = np.ones((n_particles, 3))

      # Output array rest of cases
      else:
        data = np.ones(n_particles)                                             

      del(n_particles)                                                          
      gc.collect
    
    # DarkMatter Mass case
    if itype==1 and att=='Mass':                                                
      data      *= self.m_DM                                                    
      data.dtype = [('Mass', data.dtype)]                                       
      return data

    # The rest
    else:
      """
      # There are failures when the specified itype is not present the attribute,
       so a it is used a try-except estructure
      """

      try:                                                                      
        count = 0
        # Loop over each file and extract the data.
        for i in range(self.nfiles):                                            
        
          f      = h5py.File('%s/%s.%i.hdf5'%(self.path, self.subpath, i), 'r') 
          tmp    = f['PartType%i/%s'%(itype, att)][...]                         
          data[count:count+len(tmp)] = tmp                                      
          count += len(tmp)
        
        """
        In the paper (arXiv: 1706.09899, part 4.1) repeat these calcs for each
        nfile (i.e. i value) but they didn't save in any variable, they overwrite
        all these variables, so it's fast if we do only in the last one, or 
        first or whatever.

        I choose the last one.
        """
	
        aexp = f['PartType%i/%s'%(itype, att)].attrs.get('aexp-scale-exponent')
        hexp = f['PartType%i/%s'%(itype, att)].attrs.get('h-scale-exponent')

        f.close()
        del(f, tmp, count, i)                                                  
        gc.collect()
        
        # convert comovil to physical units and eliminate h multiplications
        if att != 'ParticleIDs' and data.dtype != np.int32 and data.dtype != np.int64:                 
          data  = np.multiply(data, self.a**aexp * self.h**hexp, dtype='f8')

        del(aexp, hexp, tmp)                                                    
        gc.collect()

      # In case there are no particles of that itype in this snapshot
      except KeyError:                                                          
        if att == 'Coordinates' or att == 'Velocity':                           
          data = np.ones(3)*np.nan
        else:
          data = np.array([np.nan])
      
      finally:
        if att == 'Coordinates' or att == 'Velocity':                           
          data.dtype = [(att+'_%i'%(i), data.dtype) for i in [0, 1, 2]]

        elif att == 'ParticleIDs':                                              
          data.dtype = [('ParticleIDs', '<u8')]

        else:                                                                 
          data.dtype = [(att, data.dtype)]

        return data

  def periodicity(self, array, point, center=False):
    """
    A function that allows particles to be centred around the coordinates 
    of a point or a particle, based on the periodicity property of the universe.
    """

    # Particle as a point of reference: 
    if point.dtype == np.dtype([('ParticleIDs', '<u8'), ('Coordinates_0', '<f8'), 
                                ('Coordinates_1', '<f8'), 
                                ('Coordinates_2', '<f8'), 
                                ('Mass', '<f8'), ('itype', 'i1')]):             
      for i in [0, 1, 2]:                                                       
        i_point = point['Coordinates_%i'%(i)]
        # Particle as a reference point
        array['Coordinates_%i'%(i)] -= i_point                                  
        
        # For all particle beyond L/2...
        mask = array['Coordinates_%i'%(i)] > self.boxsize/2    
        # ... is placed on the other side                 
        array['Coordinates_%i'%(i)] -= mask.astype(np.int)*self.boxsize         
        del(mask)                                                               
        gc.collect()
        # Fot all particle beyond -L/2...
        mask = array['Coordinates_%i'%(i)] < -self.boxsize/2   
        # ... is placed on the other side                 
        mask = mask.astype(np.int)*self.boxsize
        array['Coordinates_%i'%(i)] += mask                                     
        
        # IF center = False -> the coordinates are retrieved
        if not(center):                                                         
          array['Coordinates_%i'%(i)] += i_point
        
        del(mask, i_point)                                                               
        gc.collect()

      del(i, point)
      gc.collect()
      return array

    # Point as reference point case (The rest is the same)
    else:                                                                       
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

  # Function to obtain the particle data of a certain sector
  def particles_prop(self, att=None, itype=None, gn=None, sgn=None):            
    """
      Possibility to upload a quantity of data directly. Depends on the need.
      att   --> specific attribute
      itype --> specific particle-type
      gn    --> specific GroupNumber
      sgn   --> specific SubGroupNumber
    """

    if itype == None and att == None:
      """
        For all particles: ID, coordinates, mass and itype attribute
      """
      with h5py.File('%s/%s.%i.hdf5'%(snap_0.path, snap_0.subpath, 0), 'r') as f:
        n_particles = f['Header'].attrs.get('NumPart_Total')                    
        DF          = np.ones(np.sum(n_particles), 
	      		      dtype=np.dtype([('ParticleIDs', 'u8'), ('Coordinates_0', '<f8'), 
                                              ('Coordinates_1', '<f8'), 
                                              ('Coordinates_2', '<f8'), 
                                              ('Mass', '<f8'), ('itype', 'i1')]))
	
      count = 0
      for itype, i in zip([0, 1, 4, 5], [0, 1, 2, 3]):		                      

        if n_particles[i] == 0:
          continue

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

      del(count, n_particles, att, itype, i)                                    
      gc.collect()


    elif itype == None:
      """
        For all particles: the attribute selected
      """
      dicc_dtype = {'ParticleIDs':[('ParticleIDs'  , 'u8')] , 
                    'Coordinates':[('Coordinates_0', '<f8') ,
                                   ('Coordinates_1', '<f8') , 
                                   ('Coordinates_2', '<f8')],
                    'Mass'       :[('Mass'         , '<f8')],
                    'itype'      :[('itype'        , 'i1')]}


      with h5py.File('%s/%s.%i.hdf5'%(self.path, self.subpath, 0), 'r') as f: 
        n_particles = f['Header'].attrs.get('NumPart_Total')                      
        DF          = np.ones(np.sum(n_particles), dtype=np.dtype(dicc_dtype[att]))
	
      count = 0
      for itype, i in zip([0, 1, 4, 5], [0, 1, 2, 3]):			                    

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
          DF[att][count:count+n_particles[i]] = self.read_dataset(itype= itype,att= att)

        count += n_particles[i]

      del(count, n_particles, itype, i)                                        
      gc.collect()


    else:
      """
        For the itype particles selected: the selected attribute
      """
      with h5py.File('%s/%s.%i.hdf5'%(snap_0.path, snap_0.subpath, 0), 'r') as f: 
        n_particles = f['Header'].attrs.get('NumPart_Total')[itype]            
        DF = np.ones(n_particles, dtype=np.dtype([('ParticleIDs', 'u8'),
                                                  ('Coordinates_0', '<f8'), 
                                                  ('Coordinates_1', '<f8'), 
                                                  ('Coordinates_2', '<f8'), 
                                                  ('Mass', '<f8'), ('itype', 'i1')]))   
	
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

        del(count, n_particles, att)                                            
        gc.collect()


    if gn != None and sgn != None and self.snapnum != 0:
      """  
        The data obtained is specified for a gn and sgn selected. 
        NOTE: The zeroth snapshot does not GroupNumber nor SubGroupNumber                    
      """
      gns  = self.read_dataset(itype, 'GroupNumber')['GroupNumber']
      sgns = self.read_dataset(itype, 'SubGroupNumber')['SubGroupNumber']
      mask = np.logical_and(gns == gn, sgns == sgn)                         
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

def Galaxy_to_past(GalaxyID, snap_1, snap_2):                                   
  """
    # Function to obtain center and the radio of the sphere with all 
    the particles at snap_2 of the galaxy at snap_1

    snap_1.  --> Galaxy snapshot
    snap_2.  --> Past snapshot
    GalaxyID --> ID of the galaxy at snap_1
    
    We only consider to use the darkmatter 
  """

  # GroupNumber and SubGroupNumber of the Galaxy
  mask_gn_sgn = snap_1.catalogue['GalaxyID'] == GalaxyID                        
  gn, sgn     = snap_1.catalogue[['GroupNumber', 
                                  'SubGroupNumber']][mask_gn_sgn][0]           

  del(mask_gn_sgn)                                                        
  gc.collect()
  
  gns    = snap_1.read_dataset(itype= 1, att='GroupNumber')['GroupNumber']
  sgns   = snap_1.read_dataset(itype= 1, att='SubGroupNumber')['SubGroupNumber']
  mask_1 = np.logical_and(gns == gn, sgns == sgn)

  # Array with the IDs of the galaxy at snap_1
  n_particles   = np.sum(mask_1)                    
  ParticleIDs_1 = snap_1.read_dataset(itype= 1, att='ParticleIDs')['ParticleIDs'][mask_1] 
    
  del(gns, sgns, mask_1)                                                        
  gc.collect()

  # All IDs Particles at snap_2
  ParticleIDs_2 = snap_2.read_dataset(itype = 1, att='ParticleIDs')['ParticleIDs'] 
  # Mask to looking where are the galaxy's particles on snapshot 2.
  mask_2        = np.in1d(ar1 = ParticleIDs_2, ar2 = ParticleIDs_1)             

  del(ParticleIDs_1, ParticleIDs_2)                                             
  gc.collect()
  
  # Their Coordinates
  Coordinates = snap_2.read_dataset(itype = 1, att='Coordinates')[mask_2][:,0]  
  del(mask_2)
  gc.collect()
  Coordinates = snap_2.periodicity(Coordinates, Coordinates[0].copy())            
  
  # Their mass
  Mass_T      = snap_2.m_DM*n_particles				            
  del(n_particles)
  gc.collect()
  
  # Center of mass
  Mass_center = np.array([np.sum(Coordinates['Coordinates_%i'%(i)]*snap_2.m_DM) for i in [0, 1, 2]])/Mass_T
  del(Mass_T)
  gc.collect()
  
  # the radios of the galaxy's particles with the new mass center
  Radios = np.sum([(Coordinates['Coordinates_%i'%(i)]-Mass_center[i])**2 for i in [0, 1, 2]], axis=0)  
  del(Coordinates)
  gc.collect()

  # Maximum radio of the 90$ nearest particles.
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

  with h5py.File('%s/%s.%i.hdf5'%(snap.path, snap.subpath, 0), 'r') as f:   
    n_particles = f['Header'].attrs.get('NumPart_Total')                        
    Radios_2    = np.ones(np.sum(n_particles), dtype=np.dtype('<f8'))           

  # radius of gas and DM particles with respect the center of mass.
  count = 0
  for i in [0, 1]:                                                              
    Coordinates = snap.read_dataset(itype= i, att='Coordinates')
    Coordinates = snap.periodicity(Coordinates, center)
    Radios_2[count:count+n_particles[i]] = np.sum([(Coordinates['Coordinates_%i'%(i)] \
                                                    -center[i])**2 for i in [0, 1, 2]], axis=0)[:,0]

    del(Coordinates)
    gc.collect()
      
    count += n_particles[i]
  
  del(i)
  gc.collect()

  # Total mass, density and overdensity inside of the sphere
  Mass = np.sum(snap.particles_prop(att='Mass')['Mass'][Radios_2 <= Radio])     
  rho  = Mass*1.989e43/((4/3)*np.pi*((Radio)**(1/2)*3.08e24)**3)                
  overrho = (rho - snap.rho_matter)/snap.rho_matter                             
  return overrho, rho, Mass



"""	# AngularMoment      """


def AngularMoment(snap, Radio, center):
  """
	Function to obtain the total vector angular momentum of the sphere
	given as input (radio and center) at the snapshot also selected (snap)
  """
  with h5py.File('%s/%s.%i.hdf5'%(snap.path, snap.subpath, 0), 'r') as f:   
    n_particles = f['Header'].attrs.get('NumPart_Total')                        
    Radios_2    = np.ones(np.sum(n_particles), dtype=np.dtype('<f8'))          

  count = 0
  # radius of gas and DM particles with respect the center of mass.
  for i in [0, 1]:                                                              
    Coordinates = snap.read_dataset(itype= i, att='Coordinates')		            
    Coordinates = snap.periodicity(Coordinates, center)
    Radios_2[count:count+n_particles[i]] = np.sum([(Coordinates['Coordinates_%i'%(i)] \
                                                    -center[i])**2 for i in [0, 1, 2]], axis=0)[:,0]

    del(Coordinates)
    gc.collect()
      
    count += n_particles[i]
  
  del(i)
  gc.collect()
  
  count = 0
  count_angular = 0
  angulars = np.ones((np.sum(Radios_2 <= Radio), 3), dtype=np.dtype('<f8'))
  # Calculus the vectorial product between the coordinates and the velocity of the particles inside the region
  for i in [0, 1]:                                                              
    mask = [Radios_2[count:count+n_particles[i]] <= Radio][0]
    Coordinates = snap.read_dataset(itype= i, att='Coordinates')[..., 0][mask]
    Coordinates = snap.periodicity(Coordinates, center)
    Velocity    = snap.read_dataset(itype= i, att='Velocity')[..., 0][mask]

    angulars[count_angular:count_angular + np.sum(mask)][..., 0] = Coordinates['Coordinates_1']*Velocity['Velocity_2'] \
                                                                 - Coordinates['Coordinates_2']*Velocity['Velocity_1']
    angulars[count_angular:count_angular + np.sum(mask)][..., 1] = Coordinates['Coordinates_2']*Velocity['Velocity_0'] \
                                                                 - Coordinates['Coordinates_0']*Velocity['Velocity_2']
    angulars[count_angular:count_angular + np.sum(mask)][..., 2] = Coordinates['Coordinates_0']*Velocity['Velocity_1'] \
                                                                 - Coordinates['Coordinates_1']*Velocity['Velocity_0']
   
    del(Coordinates, Velocity)
    gc.collect()
      
    count += n_particles[i]
    count_angular += np.sum(mask)

    # times DM mass for DM particles
    if i == 1:								     
      angulars[count_angular:count_angular + np.sum(mask)] *= snap.m_DM
    
    # times gas mass for gas particles
    elif i == 0:								     
      mass = snap.read_dataset(itype= i, att='Mass')['Mass'][mask][0]
      angulars[count_angular:count_angular + np.sum(mask)] *= mass
      del(mass)
    
    del(mask)

  del(i, Radios_2)
  gc.collect()
  
  # Total Angular Momentum
  angular = np.array([np.sum(angulars[i]) for i in [0, 1, 2]])
  return angular
