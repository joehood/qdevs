
"""Generic LIM Model."""



class LimNode(object):

    def __init__(self, index):
        pass


class LimDevice(object):

    def __init__(self):
        pass


class LimNodeDevice(LimDevice):

    def __init__(self):
        pass


class LimBranchDevice(LimDevice):

    def __init__(self):
        pass


class LimNetwork(object):

    def __init__(self):
        self.nodes = []

