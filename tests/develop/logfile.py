class Pymatmc2LogFile:

    def __init__(self, path):
        self.path = path

    def log(self, message): 
        print(message)