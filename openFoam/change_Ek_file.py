path = '/home/sq-wm/OpenFOAM/sq-wm-v2306/run/calc10_struct/Ek/5/'

with open(path + 'Ek.xy', 'r') as fp:
    data = fp.read().split('\n')

file = open(path + 'Ek_new.xy', "a+")
for i in range(0, len(data)):
    num = data[i].split(' ')
    nums = []
    for j in range(0, len(num)):
        if len(num[j]) > 0:
            nums.append(num[j])
    if len(nums) == 2:
        num1 = nums[0].split('.')
        num2 = nums[1].split('.')
        file.write(str(num1[0]) + ',' + str(num1[1]) + '\t' + str(num2[0]) + ',' + str(num2[1]) + '\n')
file.close()
