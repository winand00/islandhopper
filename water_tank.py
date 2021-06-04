
hydrogen_energy_density = 120 * (10 ** 6)

def water_tank(power_fc, time_to):
    """Calculates water mass in kg and volume in L"""
    hydrogen_mass_flow = power_fc / hydrogen_energy_density
    water_mass_flow = hydrogen_mass_flow * 9
    water_mass = water_mass_flow * time_to
    water_volume = (water_mass / 997) * 1000
    return water_mass, water_volume

water_mass, water_volume = water_tank(1300000, 50)
print(f"Water mass  : {water_mass} \n"
      f"Water volume: {water_volume}")