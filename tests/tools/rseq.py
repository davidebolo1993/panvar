#!/usr/bin/env python3
"""A deterministic random DNA sequence, the same on every platform.

    rseq.py <length> <seed>

The suites first drew these with awk's srand/rand, whose generator differs between awk implementations
(macOS awk uses BSD random(), GNU awk and mawk their own), so a fixture's sequence, and the positions a
check expects, changed with the machine. This reproduces the BSD random() macOS awk calls (the additive
feedback generator of degree 31, seeded by the Park-Miller recurrence), so the fixtures keep the sequences
the expected values were written for.
"""
import sys

DEG, SEP = 31, 3


def seed_step(x):
    # Park-Miller "minimal standard" step on a signed 32-bit value (Schrage's method), as BSD libc seeds it
    x &= 0xFFFFFFFF
    if x >= 2**31:
        x -= 2**32
    if x == 0:
        x = 123459876
    hi, lo = divmod(x, 127773)
    x = 16807 * lo - 2836 * hi
    if x < 0:
        x += 0x7FFFFFFF
    return x & 0xFFFFFFFF


def rseq(n, seed):
    state = [seed & 0xFFFFFFFF] + [0] * (DEG - 1)
    for i in range(1, DEG):
        state[i] = seed_step(state[i - 1])
    front, rear = SEP, 0

    def draw():
        nonlocal front, rear
        state[front] = (state[front] + state[rear]) & 0xFFFFFFFF
        value = state[front] >> 1
        front += 1
        if front >= DEG:
            front, rear = 0, rear + 1
        else:
            rear += 1
            if rear >= DEG:
                rear = 0
        return value

    for _ in range(10 * DEG):
        draw()
    # awk: substr("ACGT", int(rand() * 4) + 1, 1), rand() = random() / 2^31
    return "".join("ACGT"[draw() >> 29] for _ in range(n))


if __name__ == "__main__":
    sys.stdout.write(rseq(int(sys.argv[1]), int(sys.argv[2])))
