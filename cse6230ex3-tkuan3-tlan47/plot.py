# import numpy as np
import matplotlib.pyplot as plt
import json


file_idx = '133382'
# print('omp_data_'+file_idx+'.json')
with open('omp_data_'+file_idx+'.json') as f:
    omp_data = json.load(f)

# print(omp_data)

# print('acc_data_'+file_idx+'.json')
with open('acc_data_'+file_idx+'.json') as f:
    acc_data = json.load(f)

# print(data)

acc_wr = []
acc_n = []
omp_wr = []
omp_n = []

for ele in omp_data:
    if ele['precision']=='float':
        omp_wr.append(float(ele['charge interactions per second']))
        omp_n.append(int(ele['num charges']))
    else:
        continue
        # d_wr.append(float(ele['charge interactions per second']))
        # d_n.append(int(ele['num charges']))

for ele in acc_data:
    if ele['precision']=='float':
        acc_wr.append(float(ele['charge interactions per second']))
        acc_n.append(int(ele['num charges']))
    else:
        continue

plt.plot(omp_n, omp_wr, '.', color='orange', label='omp')
plt.plot(acc_n, acc_wr, '.', color='blue', label='acc')
plt.xlabel("number of charged particles")
plt.ylabel("work rate (charged particle interactions per second)")
plt.title('work rate-number of charged particles (float)')
plt.xscale("log")
plt.yscale("log")
plt.grid(True, which="both")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('figure_float_'+file_idx+'.png')


plt.clf()
acc_wr = []
acc_n = []
omp_wr = []
omp_n = []

for ele in omp_data:
    if ele['precision']=='double':
        omp_wr.append(float(ele['charge interactions per second']))
        omp_n.append(int(ele['num charges']))
    else:
        continue

for ele in acc_data:
    if ele['precision']=='double':
        acc_wr.append(float(ele['charge interactions per second']))
        acc_n.append(int(ele['num charges']))
    else:
        continue

plt.plot(omp_n, omp_wr, '.', color='orange', label='omp')
plt.plot(acc_n, acc_wr, '.', color='blue', label='acc')
plt.xlabel("number of charged particles")
plt.ylabel("work rate (charged particle interactions per second)")
plt.title('work rate-number of charged particles (double)')
plt.xscale("log")
plt.yscale("log")
plt.grid(True, which="both")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('figure_double_'+file_idx+'.png')


# with open("figure_two_131431.json") as f:
#     data = json.load(f)

# f_error = []
# f_syps = []
# f_steps = []
# d_error = []
# d_syps = []
# d_steps = []

# for ele in data:
#     if ele['precision']=='float':
#         f_error.append(float(ele['error']))
#         f_syps.append(float(ele['simulated years per second']))
#         f_steps.append(int(ele['n_steps']))
#     else:
#         d_error.append(float(ele['error']))
#         d_syps.append(float(ele['simulated years per second']))
#         d_steps.append(int(ele['n_steps']))

# plt.plot(f_syps, f_error, 'o-', color='orange')
# plt.plot(d_syps, d_error, '.-', color='blue')

# for i in range(len(f_steps)):
#     plt.text(f_syps[i], f_error[i], f_steps[i])
#     plt.text(d_syps[i], d_error[i], d_steps[i])

# plt.xscale("log")
# plt.yscale("log")
# plt.xlabel("simulated years per second")
# plt.ylabel("error")
# plt.grid(True, which="both")
# plt.tight_layout()
# plt.savefig('figure_two_131431.png')