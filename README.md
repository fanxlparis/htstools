# entropy

Python script to compute the **entropy** of a bytestream

---

Assumptions:
- The bytestream is given as a file.
- The value of each byte in the file is interpreted as a symbol s in [0,255].
- The symbols are assumed to be statistically independent, i.e. P(s(i)|s(i-1)) = P(s(i)), with s(i) being the symbol at position i >= 1.
