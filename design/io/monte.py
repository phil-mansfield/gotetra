from __future__ import division
from random import *

doors = ["red", "blue", "green"]
contents = ["eva", "mary", "carmen"]

# This function returns None if the trial stopped early, (True, False) if
# you would win by switching, and (False, True) if you would win by staying.
def red_mary_trial():
    # The host hides the contents
    shuffle(contents)
    # If this isn't a red mary trial, stop the trial
    if doors.index("red") != contents.index("mary"): return None
    # Choose two doors.
    my_doors = sample(doors, 2)
    my_contents = [contents[doors.index(color)] for color in my_doors]
    # If mary isn't one of the two doors you chose, stop the trial
    if "mary" not in my_contents: return None

    # The host chooses mary
    my_contents.remove("mary")

    # report the results of swithcing and staying
    stay_result = my_contents[0]
    switch_result = "eva" if stay_result == "carmen" else "carmen"

    return stay_result == "carmen", switch_result == "carmen"

if __name__ == "__main__":
    trials = 100 * 1000
    switch_wins, stay_wins = 0, 0

    for i in xrange(trials):
        res = red_mary_trial()
        if res is None: continue
        switch_res, stay_res = res
        if switch_res: switch_wins += 1
        if stay_res: stay_wins += 1

    print "Switch win rate is %.4g" % (switch_wins / (switch_wins + stay_wins))
    print "Stay win rate is   %.4g" % (stay_wins / (switch_wins + stay_wins))
