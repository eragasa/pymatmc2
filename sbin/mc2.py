#!/apps/python/3.6-conda5.2/bin/python
import click
import os
import sys
from crontab import CronTab
# thVis might be taken care of with a pip install
sys.path.append('/users/PAA0028/eragasa/repos/pymatmc2/src')
sys.path.append('/users/PAA0028/eragasa/repos/mexm-base/src')
from pymatmc2 import MultiCellMonteCarlo
from pymatmc2 import Pymatmc2Configuration
# might have to run this command first
# export LC_ALL=en_US.UTF-8 



# cron environment read .bashrc or .bash_profile
os.environ['VASP_POTPAW_GGA'] = '/users/PAA0028/eragasa/usr/local/vasp/potpaw/potpaw_PBE.54'
 

PYMATMC2_CONTINUE_CMD = "{} --continue {}".format(__file__, os.getcwd())
PYMATMC2_STOP_CMD = "python {} --stop".format(os.path.abspath(__file__))
PYMATMC2_CONTINUE_DESC = "pymatmc2_continue"

@click.command()
@click.option('--start', 'start_option', flag_value='start')
@click.option('--restart', 'start_option', flag_value='restart')
@click.option('--continue', 'start_option', flag_value='continue')
@click.option('--stop', 'start_option', flag_value='stop')
@click.option('--none', 'start_option', flag_value='none', default=True)
@click.argument('path', default=os.getcwd())
def main(start_option, path):
    start_options = {
        'start': pymatmc2_start,
        'restart': pymatmc2_restart,
        'continue': pymatmc2_continue,
        'stop': pymatmc2_stop,
        'none': pymatmc2_none
    }
    start_options[start_option](path=path)

def pymatmc2_start(path):
    kwargs_mc2 = {
        'configuration_path':'pymatmc2.config',
        'results_path':'results',
        'logfile_path':'pymatmc2.log',
        'simulations_path':'simulations',
        'is_restart':False
    }
    
    o_mc2 = MultiCellMonteCarlo(**kwargs_mc2)
    o_mc2.run()    
    
    schedule_cron_job(
        command=PYMATMC2_CONTINUE_CMD,
        description=PYMATMC2_CONTINUE_DESC
    )

def pymatmc2_restart(path):

    os.chdir(path)
    
    kwargs_mc2 = {
        'configuration_path':'pymatmc2.config',
        'results_path':'results',
        'logfile_path':'pymatmc2.log',
        'simulations_path':'simulations',
        'is_restart':True
    }
    o_mc2 = MultiCellMonteCarlo(**kwargs_mc2)
    is_max_iterations = o_mc2.run()    

    if is_max_iterations:
        msg = 'is_max_iterations:{}'
        msg = msg.format(is_max_iterations)
        print(msg)
        pymatmc2_stop(path=path)
    else:
        schedule_cron_job(
	    command=PYMATMC2_CONTINUE_CMD,
            description=PYMATMC2_CONTINUE_DESC
        )
        

def pymatmc2_continue(path):
    os.chdir(path)
    
    kwargs_mc2 = {
        'configuration_path':'pymatmc2.config',
        'results_path':'results',
        'logfile_path':'pymatmc2.log',
        'simulations_path':'simulations',
        'is_restart':True
    }
    o_mc2 = MultiCellMonteCarlo(**kwargs_mc2)
    assert isinstance(o_mc2.configuration, Pymatmc2Configuration)
    is_max_iterations = o_mc2.run()    

    if is_max_iterations:
        msg = 'is_max_iterations:{}'
        msg = msg.format(is_max_iterations)
        print(msg)
        pymatmc2_stop(path=path)

def pymatmc2_stop(path):
    remove_cron_job(
        command=PYMATMC2_CONTINUE_CMD,
        description=PYMATMC2_CONTINUE_DESC
    )

def pymatmc2_none(path):
    msg = (
        "A description of the current option for pymatmc2\n"
        "mc2.py --start\n"
        "mc2.py --continue\n"
        "mc2.py --stop\n"
    )
    print(msg)

    user = os.environ['USER']
    my_cron = CronTab(user=user)
    for job in my_cron:
        print(job)


def schedule_cron_job(command, description):
    user = os.environ['USER']
    my_cron = CronTab(user=user)

    #create the cronjob
    job = my_cron.new(
        command=PYMATMC2_CONTINUE_CMD,
        comment=PYMATMC2_CONTINUE_DESC
    )
    job.minute.every(1)

    #write the crontab file
    my_cron.write()

    print('job_is_valid:', job.is_valid())

def remove_cron_job(command, description):
    user = os.environ['USER']
    my_cron = CronTab(user=user)
    for job in my_cron:
        if job.comment == PYMATMC2_CONTINUE_DESC:
            my_cron.remove(job)
            my_cron.write()

if __name__ == "__main__":
    main()
