import numpy as np
import matplotlib.pyplot as plt

def func1(x):
    return np.sin(np.power(x, 2))

def func2(x):
    return np.cos(np.sin(x))

def func3(x):
    return np.exp(np.sin(np.cos(x)))

def func4(x):
    return np.log(x + 3)

def func5(x):
    return np.power((x + 3), (1 / 2))

def der1(x): # derivative of sin(x^2)
    return 2 * x * np.cos(np.power(x, 2))

def der2(x): # derivative of cos(sin(x))
    return -np.sin(np.sin(x)) * np.cos(x)

def der3(x): # derivative of exp(sin(cos(x)))
    return -np.exp(np.sin(np.cos(x))) * np.cos(np.cos(x)) * np.sin(x)

def der4(x): # derivative of ln(x+3)
    return 1 / (x + 3)

def der5(x): # derivative of (x+3)^0,5
    return 1 / (2 * np.power(x + 3, 0.5))

def method1(x, h, func):
    return (func(x + h) - func(x)) / h

def method2(x, h, func):
    return (func(x) - func(x - h)) / h

def method3(x, h, func):
    return (func(x + h) - func(x - h)) / (2 * h)

def method4(x, h, func):
    frac1 = 4 / 3 * (func(x + h) - func(x - h)) / (2 * h)
    frac2 = 1 / 3 * (func(x + 2 * h) - func(x - 2 * h)) / (4 * h)
    return frac1 - frac2

def method5(x, h, func):
    frac1 = 3 / 2 * (func(x + h) - func(x - h)) / (2 * h)
    frac2 = 3 / 5 * (func(x + 2 * h) - func(x - 2 * h)) / (4 * h)
    frac3 = 1 / 10 * (func(x + 3 * h) - func(x - 3 * h)) / (6 * h)
    return frac1 - frac2 + frac3

def createStepList():
    stepList = []
    for n in np.arange(1, 22, 1):
        stepList.append(2 / (np.power(2, n)))
    return stepList

def calcError(x, method, func, der, step):
    return method(x, step, func) - der(x)

def createErrorAxis(x, method, func, der, stepList):
    error_axis = []
    for i in stepList:
        error_axis.append(calcError(x, method, func, der, i))
    return error_axis
        

def createPlot(x, methodList, func, der, stepList, funcName,file):
    step_axis = stepList
    error_axis = []

    for method in methodList:
        error_axis.append(createErrorAxis(x, method, func, der, stepList))
    error_axis = np.abs(error_axis)

    plt.xlabel('step (h)')
    plt.ylabel('Error')
    plt.title(funcName)

    plt.semilogx()
    plt.semilogy()
    plt.gca().invert_xaxis()

    for i in error_axis:
        plt.plot(step_axis, i, marker='+', markersize=4)

    LegendList = ['$\\frac{f(x+h) - f(x)}{h}$',\
            
            '$\\frac{f(x) - f(x-h)}{h}$',\
            
            '$\\frac{f(x+h) - f(x-h))}{2h}$',\
            
            '$\\frac{4}{3}\\cdot\\frac{f(x+h) - f(x-h))}{2h}\
                - \\frac{1}{3}\\cdot\\frac{f(x+2h) - f(x-2h))}{4h}$',\
            
            '$\\frac{3}{2}\\cdot\\frac{f(x+h) - f(x-h))}{2h} - \\frac{3}{5}\\cdot\\frac{f(x+2h)\
                - f(x-2h))}{4h} + \\frac{1}{10}\\cdot\\frac{f(x+3h) - f(x-3h))}{6h}$'\
            ]
    plt.legend(LegendList, loc="upper right", fontsize="7")
    plt.savefig(file, dpi=800)
    plt.clf()

funcList = [func1, func2, func3, func4, func5]
derList = [der1, der2, der3, der4, der5]
methodList = [method1, method2, method3, method4, method5]
stepList = createStepList()

folder = 'graphs/'
funcName = ['$sin(x^2)$', '$cos(sin(x))$', '$exp(sin(cos(x)))$', '$ln(x+3)$', '$(x+3)^{0.5}$']
for i in range(5):
    createPlot(2, methodList, funcList[i], derList[i], stepList, funcName[i], folder+'graph'+str(i)+'.jpg')
