import click
import os
from crontab import CronTab

PYMATMC2_CONTINUE_CMD = "python {} --continue".format(os.path.abspath(__file__))
PYMATMC2_STOP_CMD = "python {} --stop".format(os.path.abspath(__file__))
PYMATMC2_CONTINUE_DESC = "pymatmc2_continue"

@click.command()
@click.option('--start', 'start_option', flag_value='start')
@click.option('--continue', 'start_option', flag_value='continue')
@click.option('--stop', 'start_option', flag_value='stop')
@click.option('--none', 'start_option', flag_value='none', default=True)
def main(start_option):
    start_options = {
        'start': pymatmc2_start,
        'continue': pymatmc2_continue,
        'stop': pymatmc2_stop,
        'none': pymatmc2_none
    }
    start_options[start_option]()

def pymatmc2_start():
    print(PYMATMC2_CONTINUE_CMD)
    schedule_cron_job(
        command=PYMATMC2_CONTINUE_CMD,
        description=PYMATMC2_CONTINUE_DESC
    )
    print('start')

def pymatmc2_continue():
    print('continue')

def pymatmc2_stop():
    remove_cron_job(
        command=PYMATMC2_CONTINUE_CMD,
        description=PYMATMC2_CONTINUE_DESC
    )
    print('stop')

def pymatmc2_none():
    msg = (
        "A description of the current option for pymatmc2\n"
        "mc2.py --start\n"
        "mc2.py --continue\n"
        "mc2.py --stop\n"
    )
    print(msg)

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

def remove_cron_job(command, description):
    user = os.environ['USER']
    my_cron = CronTab(user=user)
    for job in my_cron:
        if job.comment == PYMATMC2_CONTINUE_DESC:
            my_cron.remove(job)
            my_cron.write()

if __name__ == "__main__":
    main()