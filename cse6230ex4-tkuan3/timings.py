
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_json("timings.json", orient="records")
plt.semilogx(df["number of particles"],df["particle timesteps per second"])
plt.ylim(bottom=0)
plt.gca().set(
    title="Particle timesteps per second",
    xlabel="Number of particles",
    ylabel="Particle timesteps per second"
    )
plt.savefig("timings.png")
