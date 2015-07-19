import abc
import matplotlib.pyplot as plt
import numpy as np
import numpy.random as rand

class AbstractMetropolis(object):
    __metaclass__ = abc.ABCMeta

    def sample(self, steps, T):
        snapshots = []
        for _ in xrange(steps):
            updates = self.update_candidates()
            for update in updates:
                dE = self.dE(update)
                if dE < 0 or rand.random() < np.exp(-dE/T): self.change(update)
            snapshots.append(self.state())
        return np.array(snapshot)

    @abc.abstractmethod
    def update_candidates(self): pass
    @abc.abstractmethod
    def dE(self, update): pass
    @abc.abstractmethod
    def change(self, update): pass
    @abc.abstractmethod
    def state(self): pass

class Potential1D(AbstractMetropolis):
    def __init__(self, f, x0, dx):
        self.f, self.x, self.dx = f, x0, dx
    def update_candidates(self):
        dx = self.dx * 2 * rand.random() - self.dx
        return x0 + dx
    def dE(self, update):
        return f(self.x + update) - f(self.x)
    def change(self, update):
        return self.x += update
    def state(self):
        return self.x

if __name__ == "__main__"
