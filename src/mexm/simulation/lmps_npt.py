import os, random, sys
from collections import OrderedDict

from mexm.simulation import LammpsSimulation, NptSimulation
from mexm.io.vasp import Poscar
from mexm.structure import SimulationCell
from mexm.io.lammps import LammpsStructure

class LammpsNptSimulation(LammpsSimulation, NptSimulation):
    simulation_type = 'lammps_npt'
    is_base_class = False
    thermostat_types = ['NoseHoover']

    """ Class for LAMMPS NPT (isothermal_isobaric) simulations

    This data class defines additional attributes and methods 

    Args:
        name (str): unique id for the task name being define
        path (str): the directory where this task will create
            input and output files for LAMMPS
        simulation_path(str): path for the simulation
        structure_path(str): path for the structure
        bulk_structure_name(str): name for the structure.
        temperature_initial (int): in degrees kelvin
        temperature_final (int): in degrees kelvin
        temperature_damp (int): temperature dampening, no units
        pressure_initial (int): in bars
        pressure_final (int): in bars
        pressure_damp (int): pressure dampening, no units
        drag_coefficient (float): drag coefficient, no units
        time_total (int): in picoseconds
        total_step (int): in femtoseconds
        pressure_damp (float): no units
    Attributes:
        lammps_npt_ramp_path (str)
        lammps_npt_path (str)
    """
    def __init__(
        self,
        name,
        simulation_path,
        structure_path,
        bulk_structure_name = None,
        thermostat_type = 'NoseHoover',
        temperature_initial = 0.0,
        temperature_final = 100.0,
        temperature_damp = 'Auto',
        pressure_initial = 0.0,
        pressure_final = 1000.0,
        pressure_damp = 'Auto',
        drag_coefficient = 'Auto',
        time_ramp = 10000,
        time_run = 10000,
        time_step = 3
    ):

        LammpsSimulation.__init__(self,
                                  name=name,
                                  path=simulation_path,
                                  structure_path=structure_path,
                                  bulk_structure_name=None)

        self.lammps_npt_ramp_path = os.path.join(self.path, 'lmps_npt_ramp.out')
        self.lammps_npt_path = os.path.join(self.path, 'lmps_npt.out')

        self._thermostat_type = thermostat_type
        self._temperature_initial = temperature_initial
        self._temperature_final = temperature_final
        self._temperature_damp = temperature_damp
        self._pressure_initial = pressure_initial
        self._pressure_final = pressure_final
        self._pressure_damp = pressure_damp
        self._drag_coefficient = drag_coefficient
        self._time_ramp = time_ramp
        self._time_run = time_run
        self._time_step = time_step

    @property
    def thermostat_type(self):
        return self._thermostat_type
    
    @thermostat_type.setter
    def thermostat_type(self, thermostat_type):
        if thermostat_type not in self.thermostat_types:
            msg = "{} is not a supported thermostat".format(thermostat_type)
            raise ValueError(msg)
        self._thermostat_type = thermostat_type

    @property
    def temperature_initial(self):
        return self._temperature_initial

    @temperature_initial.setter
    def temperature_initial(self, temperature):
        temperature_ = float(temperature)
        if temperature_ < 0.:
            msg = "temperature cannot be less than zero"
            raise ValueError(msg)
        self._temperature_initial = temperature_
    
    @property
    def temperature_final(self):
        return self._temperature_final
        
    @temperature_final.setter
    def temperature_final(self, temperature):
        temperature_ = float(temperature)
        if temperature_ < 0.:
            msg = "temperature cannot be less than zero"
            raise ValueError(msg)
        self._temperature_final = temperature_

    @property
    def temperature_damp(self):
        if self._temperature_damp == 'Auto':
            return self.get_recommended_temperature_damp()
        else:
            return self._temperature_damp

    @temperature_damp.setter
    def temperature_damp(self, t_damp):
        self._temperature_damp = t_damp

    @property
    def pressure_initial(self):
        return self._pressure_initial

    @pressure_initial.setter
    def pressure_initial(self, pressure):
        pressure_ = float(pressure)
        if pressure_ < 0:
            msg = "pressure cannot be negative"
            raise ValueError(msg)
        self._pressure_initial = pressure_

    @property
    def pressure_final(self):
        return self._pressure_final
    
    @pressure_final.setter
    def pressure_final(self, pressure):
        pressure_ = float(pressure)
        if pressure_ < 0:
            msg = "pressure cannot be negative"
            raise ValueError(msg)
        self._pressure_final = pressure_

    @property
    def pressure_damp(self):
        if self._pressure_damp == 'Auto':
            return self.get_recommended_pressure_damp()
        else:
            return self._pressure_damp

    @pressure_damp.setter
    def pressure_damp(self, p_damp):
        if p_damp == 'Auto':
            self._pressure_damp = 'Auto'
        else:
            self._pressure_damp = float(p_damp)

    @property
    def drag_coefficient(self):
        if self._drag_coefficient == 'Auto':
            return self._determine_drag_coefficient()
        else:
            return self._drag_coefficient

    @drag_coefficient.setter
    def drag_coefficient(self, drag_c):
        if drag_c == 'Auto':
            self._drag_coefficient = 'Auto'
        else:
            self._drag_coefficent = float(drag_c)

    @property
    def time_ramp(self):
        return self._time_ramp

    @time_ramp.setter
    def time_ramp(self, t_ramp):
        time_ramp = int(t_ramp)
        if time_ramp < 0:
            msg = "time_ramp must be greater than zero"
            raise ValueError(msg)
        self._time_ramp = time_ramp

    @property
    def time_run(self):
        return self._time_run

    @time_run.setter
    def time_run(self, t_run):
        time_run = int(t_run)
        if time_run < 0:
            msg = "time_run must be greater than zero"
            raise ValueError(msg)
        self._tim_rune = time_run

    @property
    def time_step(self):
        return self._time_step

    @time_step.setter
    def time_step(self, t_step):
        time_step = int(t_step)
        if time_step < 0:
            msg = "time step must be greater than zero"
        self._time_step = time_step

    @property
    def npt_configuration(self):
        configuration = {
            'thermostat_type':self.thermostat_type,
            'temperature_initial':self.temperature_initial,
            'temperature_final':self.temperature_final,
            'temperature_damp':self.temperature_damp,
            'pressure_initial':self.pressure_initial,
            'pressure_final':self.pressure_final,
            'pressure_damp':self.pressure_damp,
            'drag_coefficient':self.drag_coefficient,
            'time_ramp':self.time_ramp,
            'time_run':self.time_run,
            'time_step':self.time_step
        }
        return configuration


    def modify_structure(self, sc=None):
        if sc is not None:
            self.supercell=sc

    def set_npt_thermostat(self,
                           temperature=300,
                           pressure=0,
                           time_total=None,
                           time_step=None,
                           pressure_damp=None,
                           temperature_damp=None,
                           supercell=[1,1,1]):
        """
        Arguments:
            temperature (int): in degrees kelvin
            pressure (int): in bars
            time_total (int): in picoseconds
            total_step (int): in femtoseconds
            pressure_damp (float): no units
            supercell (list of int): supercell in sc1, sc2, sc3.
        """

        if not isinstance(time_total, int):
            raise TypeError('time_total must be an integer greater than zero.')
        if time_total <= 0:
            raise TypeError('time_total must be an integer greater than zero.')
        self.time_total = time_total

        if not isinstance(time_step, int):
            raise TypeError('time_step must be an integer greater than zero.')
        if not isinstance(time_step, int):
            raise TypeError('time_step must be an integer greater than zero.')
        self.time_step = time_step

        self.npt_temperature = temperature
        self.npt_pressure = pressure

        if pressure_damp is not None:
            self.npt_pressure_damp = pressure_damp
        else:
            self.npt_pressure_damp = self.get_recommended_pressure_damp()

        if temperature_damp is not None:
            self.npt_temperature_damp = temperature_damp
        else:
            self.npt_temperature_damp = self.get_recommended_temperature_damp()

        self.time_total = None
        self.time_step = None

    def get_recommended_pressure_damp(self):
        return self._determine_pressure_dampening(dt=self.time_step)

    def get_recommended_temperature_damp(self):
        return self._determine_temperature_dampening(dt=self.time_step)

    def get_task_name(self, 
                      structure, 
                      temperature,
                      pressure):
        """ get the task name

        Given a structure and a temperature, this method provides a consistent
        naming convention for a simulation name.

        Arguments:
            structure(str): structure name
            temperature(float,int): temperature in Kelvin
            pressure(float, int): pressure in GPa
        """
        if not isinstance(structure, str):
            raise TypeError('structure must be a string.')

        try:
            T = int(temperature)
            T = str(T)
        except ValueError:
            msg = 'could not cast temperature into string: {}'.format(temperature)
            raise ValueError(msg)

        try:
            P = int(pressure)
            P = str(P)
        except ValueError:
            msg = 'could not cast pressure into string: {}'.format(pressure)
            raise ValueError(msg)

        task_name = '{s}.lmps_npt_{T}K_{P}bar'.format(s=structure, T=T, P=P)
        return task_name

    def on_post(self,configuration=None):
        self._get_lattice_parameter_from_lattice_out_file(
            filename=os.path.join(self.task_directory,self.lattice_fn),
            supercell=self.supercell)

        LammpsSimulation.on_post(self,configuration=configuration)

    def write_lammps_structure_filename(self):
        _supercell = self.supercell
        _structure_filename_src \
                = self.configuration['structure']['structure_filename']
        _structure = vasp.Poscar()
        _structure.read(_structure_filename_src)

        self.structure = crystal.make_super_cell(
                structure=_structure,
                sc=_supercell)

        self.lammps_structure = LammpsStructure(obj=self.structure)

        _str_out = self.lammps_input_file_to_string()

        _lammps_filename = os.path.join(self.task_directory, filename)
        with open(_lammps_filename,'w') as f:
            f.write(_str_out)

    def lammps_input_file_to_string(self):
        str_out = "".join([\
                self._lammps_input_initialization_section(),
                self._lammps_input_create_atoms(),
                self._lammps_input_define_potential(),
                self._lammps_input_run_minimization(),
                self._lammps_input_npt_thermostat(
                    time_total = self.time_total,
                    time_step = self.time_step,
                    temperature = self.npt_temperature,
                    pressure=self.npt_pressure),
                self._lammps_input_out_section(),
                ])
        return(str_out)

    def _lammps_input_create_atoms(self):
        str_out = LammpsSimulation._lammps_input_create_atoms(self)
        str_out += "replicate {} {} {}\n".format(
                self.supercell[0],
                self.supercell[1],
                self.supercell[2])
        str_out += "change_box all x scale 1 y scale 1 z scale 1 remap\n"

        return str_out

    def _lammps_input_run_minimization(self):
        str_out = (
            '# ---- define settings\n'
            'compute eng all pe/atom\n'
            'compute eatoms all reduce sum c_eng\n'
            '# ---- run minimization\n'
            'reset_timestep 0\n'
            'fix 1 all box/relax iso 0.0 vmax 0.001\n'
            'thermo 10\n'
            'thermo_style custom step pe lx ly lz xy xz yz press pxx pyy pzz pxy pxz pyz c_eatoms\n'
            # 'thermo_style custom step pe Lx ly lz press pxx pyy pzz c_eatoms\n'
            'min_style cg\n'
            'minimize 1e-25 1e-25 5000 10000\n'
            'unfix 1\n'
            )
        return str_out

    def _lammps_input_out_section(self):
        str_out = (
            '# ---- define output variables ----\n'
            'variable natoms equal "count(all)"\n'
            'variable tot_energy equal "c_eatoms"\n'
            'variable a11 equal "xhi-xlo"\n'
            'variable a22 equal "yhi-ylo"\n'
            'variable a33 equal "zhi-zlo"\n'
            'variable tilt_xy equal "xy"\n'
            'variable tilt_xz equal "xz"\n'
            'variable tilt_yz equal "yz"\n'
            'variable tot_press equal "press"\n'
            'variable press_xx equal "pxx"\n'
            'variable press_yy equal "pyy"\n'
            'variable press_zz equal "pzz"\n'
            'variable press_xy equal "pxy"\n'
            'variable press_xz equal "pxz"\n'
            'variable press_yz equal "pyz"\n'
            '\n'
            '# ---- output ----\n'
            'print \"pypospack:output_section:begin\"\n'
            'print \"num_atoms = ${natoms}"\n'
            'print \"a11 = ${a11}\"\n'
            'print \"a22 = ${a22}\"\n'
            'print \"a33 = ${a33}\"\n'
            'print \"a12 = ${tilt_xy}\"\n'
            'print \"a13 = ${tilt_xz}\"\n'
            'print \"a23 = ${tilt_yz}\"\n'
            'print \"tot_press = ${tot_press}\"\n'
            'print \"pxx = ${press_xx}\"\n'
            'print \"pyy = ${press_yy}\"\n'
            'print \"pzz = ${press_zz}\"\n'
            'print \"pxy = ${press_xy}\"\n'
            'print \"pxz = ${press_xz}\"\n'
            'print \"pyz = ${press_yz}\"\n'
            'print \"pypospack:output_section:done\"\n'
            'print \"pypospack:lammps_sim:done\"\n'
                  )
        return str_out

 
    def _determine_temperature_dampening(self, dt):
        # A Nose-Hoover thermostat will not work well for arbitrary values of
        # Tdamp. If Tdamp is too small, the temperature can fluctuate wildly;
        # if it is too large, the temperature will take a very long time to
        # equilibrate. A good choice for many models is a Tdamp of around 100
        # timesteps.
        # Ref: http://lammps.sandia.gov/doc/fix_nh.html#fix-npt-command

        _tempdamp = dt*100
        return _tempdamp

    def _determine_pressure_dampening(self,dt):
        # A Nose-Hoover barostat will not work well for arbitrary values of
        # Pdamp. If Pdamp is too small, the pressure and volume can fluctuate
        # wildly; if it is too large, the pressure will take a very long time
        # to equilibrate. A good choice for many models is a Pdamp of around
        # 1000 timesteps.
        # Ref: http://lammps.sandia.gov/doc/fix_nh.html#fix-npt-command

        _pressdamp = dt*1000
        return _pressdamp

    def _determine_drag_coefficient(self):
        # In some cases (e.g. for solids) the pressure (volume) and/or
        # temperature of the system can oscillate undesirably when a Nose/Hoover
        # barostat and thermostat is applied. The optional drag keyword will
        # damp these oscillations, although it alters the Nose/Hoover equations.
        # A value of 0.0 (no drag) leaves the Nose/Hoover formalism unchanged.
        # A non-zero value adds a drag term; the larger the value specified,
        # the greater the damping effect. Performing a short run and monitoring
        # the pressure and temperature is the best way to determine if the drag
        # term is working. Typically a value between 0.2 to 2.0 is sufficient
        # to damp oscillations after a few periods. Note that use of the drag
        # keyword will interfere with energy conservation and will also change
        # the distribution of positions and velocities so that they do not
        # correspond to the nominal NVT, NPT, or NPH ensembles.

        _drag = 1.0
        return _drag

    def _get_random_seed(self):

        _max_size = 10000000
        _seed = random.randrange(_max_size)

        return _seed

    def _lammps_input_npt_thermostat(self,
            time_step=0.001,
            time_total=100,
            temperature=100,
            pressure=0,
            seed1=None,
            seed2=None):
        """

        Args:
            time_step (float): time increment in picoseconds
            time_total (float): time simulation time in picoseconds
            temperature (int): simulation temperature in Kelvin
            pressure (int): simulation pressure in KBar
        """

        _dt = time_step
        _time1 =time_total
        _n_time_steps = int(_time1 / _dt)
        _temp0 = temperature
        _temp1 = temperature
        _press0 = pressure
        _press1 = pressure

        _tempdamp = self._determine_temperature_dampening(dt=_dt)
        _pressdamp = self._determine_pressure_dampening(dt=_dt)
        _drag = self._determine_drag_coefficient()

        if seed1 is None:
            _seed1 = self._get_random_seed()
        else:
            _seed1 = seed1

        if seed2 is None:
            _seed2 = self._get_random_seed()
        else:
            _seed2 = seed2

        str_out = "\n".join([
            "#------------------------------------------------------------------------------",
            "# RUN THERMOSTAT",
            "# running using an NPT Nose-Hoover style thermostat",
            "#------------------------------------------------------------------------------",
            "variable tempdamp equal {tempdamp}".format(tempdamp=_tempdamp),
            "variable pressdamp equal {pressdamp}".format(pressdamp=_pressdamp),
            "",
            "timestep {dt}".format(dt=_dt),
            "# set thermo -----------------------------------------------------------------",
            "thermo 100",
            "thermo_style custom step temp pe ke etotal press lx ly lz press pxx pyy pzz pxy pxz pyz vol",
            "thermo_modify flush yes",
            "# set averaging ---------------------------------------------------------------",
            "variable boxx equal lx",
            "variable boxy equal ly",
            "variable boxz equal lz",
            "variable boxp equal press",
            "variable boxt equal temp",
            "# calculate averages every 10000 steps",
            "# set initial velocities ------------------------------------------------------",
            "reset_timestep 0",
            "velocity all create {temp} {seed} mom yes dist gaussian loop all".format(
                temp=_temp0,
                seed=_seed1),
            "# ramping the temperature",
            "# fix for Nose-Hoover style thermostat ----------------------------------------",
            "fix npt1 all npt temp {temp0} {temp1} {tempdamp} aniso 0.0 0.0 {pressdamp} drag {drag} couple xyz".format(
                temp0=_temp0,temp1=_temp1,tempdamp=_tempdamp,
                press0=_press0,press1=_press1,pressdamp=_pressdamp,
                drag=_drag),
            "fix npt1out all ave/time 1 500 500 v_boxx v_boxy v_boxz v_boxp v_boxt file lattice1.out",
            "run {n_time_steps}".format(n_time_steps=_n_time_steps),
            "unfix npt1",
            "unfix npt1out",
            "# fix for Nose-Hoover style thermostat ----------------------------------------",
            "fix npt2 all npt temp {temp0} {temp1} {tempdamp} aniso 0.0 0.0 {pressdamp} drag {drag} couple xyz".format(
                temp0=_temp1,temp1=_temp1,tempdamp=_tempdamp,
                press0=_press0,press1=_press1,pressdamp=_pressdamp,
                drag=_drag),
            "fix npt2out all ave/time 1 500 500 v_boxx v_boxy v_boxz v_boxp v_boxt file lattice2.out",
            "compute rdf2 all rdf 50",
            "fix rdf2out all ave/time 100 1 100 c_rdf2[*] file rdf2.out mode vector",
            "run {n_time_steps}\n".format(n_time_steps=_n_time_steps)
        ])

        return str_out

    def _get_lattice_parameter_from_lattice_out_file(self,filename,supercell):
        with open(filename) as f:
            lines = f.readlines()

        names = [v.strip() for v in lines[1].strip().split(" ")[1:]]
        values = [float(v) for v in lines[len(lines)-1].strip().split(" ")]

        if self.results is None:
            self.results = OrderedDict()

        _task_name = self.task_name
        _a1 = values[names.index('v_boxx')]/supercell[0]
        _a2 = values[names.index('v_boxy')]/supercell[1]
        _a3 = values[names.index('v_boxz')]/supercell[2]

        self.results['{}.{}'.format(_task_name,'a1')] = _a1
        self.results['{}.{}'.format(_task_name,'a2')] = _a2
        self.results['{}.{}'.format(_task_name,'a3')] = _a3

        return self.results
