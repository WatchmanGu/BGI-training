# -*- coding: utf-8 -*-
import time
import calendar

localtime = time.localtime(time.time())
print("本地时间为：",localtime)

localtime = time.asctime(time.localtime(time.time()))
print("本地时间：",localtime)

# 格式化成2016-03-20 11:45:39形式
print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())) 

# 格式化成Sat Mar 28 22:24:24 2018形式
print(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime())) 
  
# 将格式字符串转换为时间戳
a = "Sat Mar 28 22:24:24 2018"
print(time.mktime(time.strptime(a,"%a %b %d %H:%M:%S %Y")))

#输出某月的日历
cal = calendar.month(2018, 4)
print("以下输出2018年4月份的日历:")
print(cal)

#计算运行时间
start = time.clock()

end = time.clock()
print('Running time:%s seconds'%(end - start))