import numpy as np
from matplotlib import pyplot as plt

# MC = manufacturing cost
# VC = variable cost

#Parameters
RDTE = 64571000 / 1000000
ref_ac_MC = 7800000 * 0.5
fuel_tank_MC = 101400
battery_MC = 21697
fuel_cell_MC = 2238656
VC = (ref_ac_MC + fuel_tank_MC + battery_MC + fuel_cell_MC) / 1000000
revenue_per_ac = 10000000 / 1000000
plot_years = 3
sales_per_year = 23

#Calculate Break-Even Units (BEU) and Break-Even Revenue (BER)
BEU = 0
BER = 0
for i in range((plot_years * sales_per_year) + 1):
    if (revenue_per_ac * i) > (RDTE + VC * i):
        BEU = i
        BER = revenue_per_ac * i
        break

#Plot
fig, ax1 = plt.subplots()
ax1.plot(np.arange(0, (plot_years * sales_per_year) + 1), revenue_per_ac * np.arange(0, (plot_years * sales_per_year) + 1), label='Total Revenue')
ax1.plot(np.arange(0, (plot_years * sales_per_year) + 1), RDTE + VC * np.arange(0, (plot_years * sales_per_year) + 1), label='Total Cost')
plt.annotate('Break-Even Point', xy=(BEU + 3, BER - 30))
# ax1.plot(BEU, BER, 'kx') # label='Break-Even Point'
# ax1.plot(np.arange(0, (plot_years * sales_per_year) + 1), RDTE * np.ones((plot_years * sales_per_year) + 1), 'k--', alpha=0.5, label='Fixed Cost')

#Add BEP lines
ax1.plot([0, BEU], [BER, BER], 'k--', alpha=0.5)
ax1.plot([BEU, BEU], [0, BER], 'k--', alpha=0.5)

#Customize graph
ax2 = ax1.twiny()
ax2.set_xticks(np.arange(0, plot_years) + 1)

plt.title('Break-Even Point', weight='bold')
ax1.legend(loc='best')
ax1.grid(b=True, which='major', axis='y')
ax1.set_xlabel('Aircraft sold [-]')
ax2.set_xlabel('Years [-]')
ax1.set_ylabel('million [$]')
plt.show()

# #Plot
# fig, ax1 = plt.subplots()
# ax1.plot(np.arange(0, (plot_years * sales_per_year) + 1), (revenue_per_ac - VC) * np.arange(0, (plot_years * sales_per_year) + 1), label='Profit')
# ax1.plot(np.arange(0, plot_years * 23), RDTE * np.ones(plot_years * 23), label='RDTE')
#
# #Customize graph
# ax2 = ax1.twiny()
# ax2.set_xticks(np.arange(0, plot_years) + 1)
#
# plt.title('Break-Even Point', weight='bold')
# ax1.legend(loc='best')
# ax1.grid(b=True, which='major', axis='y')
# ax1.set_xlabel('Aircraft sold [-]')
# ax2.set_xlabel('Years [-]')
# ax1.set_ylabel('million [$]')
# plt.show()