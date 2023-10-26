import random
import math


def prime_test(N, k):
    # This is main function, that is connected to the Test button. You don't need to touch it.
    return fermat(N, k), miller_rabin(N, k)


# The mod_exp function has a time complexity O(n^3)
# The function does O(n) recursive calls with y being halved
# And each recursive call take O(n^2) because of the multiplication
# The space complexity of the algorithm is O(n) because each recursive call creates a new x, y, and N to be stored
def mod_exp(x, y, N):
    if y == 0:
        return 1
    z = mod_exp(x, math.floor(y / 2), N)  # O(n) times called recursively
    if y % 2 == 0:
        z *= z  # O(n^2)
        if z > N:
            while z > N:
                z -= N
        return z
    else:
        z = z * z * x  # O(n^2)
        if z > N:
            while z > N:
                z -= N
        return z


# The fprobability function has a time complexity of O(n^2)
# Since you are multiplying and dividing
# The space complexity is O(1)
def fprobability(k):
    k = 1 - (1 / 2 ** k)
    return k


# The mprobability function has a time complexity of O(n^2)
# Since you are multiplying and dividing
# The space complexity is O(1)
def mprobability(k):
    k = 1 - (1 / 4 ** k)
    return k


# The fermat funciton has a time complexity of O(kn^3)
# The function loops through the for loop k times
# And in each loop it calls the mod_exp function with a time complexity O(n^3)
# The space complexity is O(n) because recursion takes place in the mod_exp function
def fermat(N, k):
    for a in range(k):  # O(k)
        a = random.randint(1, N - 1)
        if mod_exp(a, N - 1, N) != 1:  # O(n^3)
            return 'composite'
    return 'prime'


# The miller_rabin function has a time complexity of O(kn^4)
# The function loops through the for loop k times
# For each k time through we run through the while loop n times
# And call the mod_exp function O(n^3)
# The space complexity is O(n) because recursion takes place in the mod_exp function
def miller_rabin(N, k):
    for a in range(k):  # O(k)
        t = N - 1
        a = random.randint(1, t)
        if mod_exp(a, t, N) != 1:  # O(n^3)
            return 'composite'
        while t % 2 == 0:  # O(n)
            t /= 2
            i = mod_exp(a, t, N)  # O(n^3)
            if i == N - 1:
                break
            else:
                if i != 1:
                    return 'composite'
    return 'prime'
