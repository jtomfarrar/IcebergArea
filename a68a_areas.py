import matplotlib.pyplot as plt
import pandas as pd

#other literature data
knowndata = [5800, 3900,2700]
knowndates = ['2020-01-01', '2020-12-15','2021-01-01']

#raw data
data2020 = {'2020-11-02': 4051.092405713962, '2020-01-01': 'Too cloudy.', '2020-01-02': 'Too cloudy.', '2020-01-03': 'Too cloudy.', '2020-01-04': 'Too cloudy.', '2020-01-05': 'Too cloudy.', '2020-01-06': 'Too cloudy.', '2020-01-07': 'Too cloudy.', '2020-01-09': 'Too cloudy.', '2020-01-08': 'Too cloudy.', '2020-01-10': 6213.117227725485, '2020-01-11': 5857.754928877559, '2020-01-13': 4628.238586071287, '2020-01-12': 'Too cloudy.', '2020-01-14': 'Too cloudy.', '2020-01-15': 5605.867054447171, '2020-01-16': 'Too cloudy.', '2020-01-17': 'Too cloudy.', '2020-01-18': 'Too cloudy.', '2020-01-19': 'Too cloudy.', '2020-01-20': 'Too cloudy.', '2020-01-22': 'Too cloudy.', '2020-01-21': 'Too cloudy.', '2020-01-23': 'Too cloudy.', '2020-01-25': 'Too cloudy.', '2020-01-24': 'Too cloudy.', '2020-01-26': 'Too cloudy.', '2020-01-27': 'Too cloudy.', '2020-01-28': 'Too cloudy.', '2020-01-29': 'Too cloudy.', '2020-01-30': 'Too cloudy.', '2020-01-31': 'Too cloudy.', '2020-02-01': 'Too cloudy.', '2020-02-02': 'Too cloudy.', '2020-02-03': 8443.293166912046, '2020-02-04': 6878.027756624643, '2020-02-05': 'Too cloudy.', '2020-02-06': 'Too cloudy.', '2020-02-07': 5718.880571688972, '2020-02-08': 5664.083355619599, '2020-02-09': 5256.874613851924, '2020-02-10': 5542.294818951318, '2020-02-11': 5031.580716256203, '2020-02-12': 'Too cloudy.', '2020-02-13': 5340.266438392604, '2020-02-14': 5073.032614178612, '2020-02-15': 'Too cloudy.', '2020-02-16': 4485.659492208006, '2020-02-17': 'Too cloudy.', '2020-02-18': 'Too cloudy.', '2020-02-19': 'Too cloudy.', '2020-02-20': 5014.021264194056, '2020-02-22': 2880.145517743842, '2020-02-21': 'Too cloudy.', '2020-02-23': 4536.63264338197, '2020-02-24': 'Too cloudy.', '2020-02-25': 5356.204493732771, '2020-02-26': 5282.548689172214, '2020-02-27': 'Too cloudy.', '2020-02-28': 'Too cloudy.', '2020-02-29': 5106.182817096922, '2020-03-02': 5282.33191833189, '2020-03-01': 'Too cloudy.', '2020-03-03': 5162.851398104309, '2020-03-04': 4628.784667942137, '2020-03-05': 3543.091451612174, '2020-03-06': 5875.17773444151, '2020-03-07': 5363.534549249233, '2020-03-08': 5066.6922713964805, '2020-03-09': 5384.949914276857, '2020-03-10': 4898.1145528506995, '2020-03-11': 'Too cloudy.', '2020-03-12': 3900.22929823568, '2020-03-13': 4782.511364219569, '2020-03-14': 2868.7132501249635, '2020-03-15': 4710.372559376201, '2020-03-16': 2885.113645932712, '2020-03-17': 2736.1122889909684, '2020-03-18': 'Too cloudy.', '2020-03-19': 5021.161168618658, '2020-03-20': 'Too cloudy.', '2020-03-21': 4847.292105718524, '2020-03-22': 5248.619220136494, '2020-03-23': 'Too cloudy.', '2020-03-24': 'Too cloudy.', '2020-03-25': 5046.491160394505, '2020-03-26': 'Too cloudy.', '2020-03-27': 4197.455603498768, '2020-03-28': 3995.304979457859, '2020-03-29': 5325.483523260838, '2020-03-30': 2619.1258417273493, '2020-03-31': 'Too cloudy.', '2020-04-01': 2256.861595944691, '2020-04-02': 8335.571795229687, '2020-04-03': 'Too cloudy.', '2020-04-05': 'Too cloudy.', '2020-04-04': 'Too cloudy.', '2020-04-06': 3620.1405400743993, '2020-04-07': 4827.870621414282, '2020-04-08': 3661.2575957767044, '2020-04-09': 5377.4140393617, '2020-04-10': 'Too cloudy.', '2020-04-11': 4571.870723467035, '2020-04-12': 7840.842866285632, '2020-04-13': 334.2482461879563, '2020-04-15': 3534.1628220110792, '2020-04-14': 3539.3654936266107, '2020-04-16': 2352.342157297852, '2020-04-18': 4910.5732179009265, '2020-04-17': 'Too cloudy.', '2020-04-19': 2928.190862669781, '2020-04-20': 759.5447573470127, '2020-04-23': 4733.162959020055, '2020-04-22': 1533.393552996588, '2020-04-21': 'Too cloudy.', '2020-04-24': 1192.8344669728765, '2020-04-25': 3256.5961129492366, '2020-04-26': 7206.858287447322, '2020-04-27': 1004.457983770016, '2020-04-28': 5340.442433902252, '2020-04-29': 3892.9788054587093, '2020-04-30': 'Too cloudy.', '2020-05-01': 2588.468438702106, '2020-05-02': 4518.489451069246, '2020-05-03': 4805.624763738858, '2020-05-04': 2346.8120583940367, '2020-05-05': 5256.508005229335, '2020-05-09': 'Too cloudy.', '2020-05-06': 1270.5888071790396, '2020-05-07': 'Too cloudy.', '2020-05-08': 3940.3782496164195, '2020-05-10': 'Too cloudy.', '2020-05-11': 2440.0978597525477, '2020-05-12': 'Too cloudy.', '2020-05-13': 6589.127230240928, '2020-05-14': 'Too cloudy.', '2020-05-15': 'Too cloudy.', '2020-05-17': 'Too cloudy.', '2020-05-16': 'Too cloudy.', '2020-05-18': 'Too cloudy.', '2020-05-19': 'Too cloudy.', '2020-05-20': 'Too cloudy.', '2020-05-21': 'Too cloudy.', '2020-05-22': 5885.74573578561, '2020-05-23': 'Too cloudy.', '2020-05-24': 'Too cloudy.', '2020-05-25': 'Too cloudy.', '2020-05-26': 'Too cloudy.', '2020-05-27': 'Too cloudy.', '2020-05-28': 'Too cloudy.', '2020-05-29': 1050.312833611482, '2020-05-30': 'Too cloudy.', '2020-05-31': 2353.439913886321, '2020-06-01': 'Too cloudy.', '2020-06-02': 'Too cloudy.', '2020-06-03': 2993.172228352871, '2020-06-04': 'Too cloudy.', '2020-06-05': 1506.7864832335908, '2020-06-06': 'Too cloudy.', '2020-06-07': 'Too cloudy.', '2020-06-08': 'Too cloudy.', '2020-06-09': 'Too cloudy.', '2020-06-10': 'Too cloudy.', '2020-06-11': 'Too cloudy.', '2020-06-12': 'Too cloudy.', '2020-06-13': 'Too cloudy.', '2020-06-14': 'Too cloudy.', '2020-06-15': 5123.097157907451, '2020-06-16': 'Too cloudy.', '2020-06-17': 4408.071402944288, '2020-06-18': 'Too cloudy.', '2020-06-20': 'Too cloudy.', '2020-06-19': 2044.8450273844242, '2020-06-21': 'Too cloudy.', '2020-06-22': 'Too cloudy.', '2020-06-23': 'Too cloudy.', '2020-06-24': 5074.976825313974, '2020-06-25': 'Too cloudy.', '2020-06-26': 1863.9918244136527, '2020-06-27': 'Too cloudy.', '2020-06-28': 'Too cloudy.', '2020-06-29': 'Too cloudy.', '2020-06-30': 'Too cloudy.', '2020-07-01': 'Too cloudy.', '2020-07-03': 'Too cloudy.', '2020-07-02': 6370.453900106993, '2020-07-06': 'Too cloudy.', '2020-07-05': 7398.443343023933, '2020-07-04': 'Too cloudy.', '2020-07-07': 1154.5632506061656, '2020-07-08': 'Too cloudy.', '2020-07-09': 'Too cloudy.', '2020-07-10': 'Too cloudy.', '2020-07-11': 'Too cloudy.', '2020-07-12': 'Too cloudy.', '2020-07-13': 'Too cloudy.', '2020-07-15': 'Too cloudy.', '2020-07-14': 'Too cloudy.', '2020-07-16': 'Too cloudy.', '2020-07-17': 'Too cloudy.', '2020-07-18': 814.2103608824009, '2020-07-19': 'Too cloudy.', '2020-07-20': 'Too cloudy.', '2020-07-21': 'Too cloudy.', '2020-07-22': 'Too cloudy.', '2020-07-23': 'Too cloudy.', '2020-07-24': 'Too cloudy.', '2020-07-25': 'Too cloudy.', '2020-07-26': 'Too cloudy.', '2020-07-27': 'Too cloudy.', '2020-07-28': 'Too cloudy.', '2020-07-29': 'Too cloudy.', '2020-07-30': 3573.2892389296917, '2020-07-31': 'Too cloudy.', '2020-08-01': 8103.879062601489, '2020-08-02': 'Too cloudy.', '2020-08-03': 735.0614252736117, '2020-08-04': 3480.8703558201937, '2020-08-05': 4198.199611654998, '2020-08-06': 2939.8869986602404, '2020-08-07': 2297.5860984252504, '2020-08-08': 'Too cloudy.', '2020-08-09': 4634.8251266090065, '2020-08-10': 485.1485579159979, '2020-08-11': 4456.557629621821, '2020-08-12': 3680.836969766394, '2020-08-13': 'Too cloudy.', '2020-08-14': 4552.704748639671, '2020-08-15': 4281.910263959917, '2020-08-16': 4636.388231612264, '2020-08-17': 876.7315284768425, '2020-08-18': 4666.599068705713, '2020-08-19': 8027.130865699742, '2020-08-20': 4302.229207854116, '2020-08-21': 'Too cloudy.', '2020-08-22': 'Too cloudy.', '2020-08-23': 'Too cloudy.', '2020-08-24': 4883.642923844344, '2020-08-25': 7737.919347364539, '2020-08-26': 4069.546987943734, '2020-08-27': 'Too cloudy.', '2020-08-28': 1275.5209230600751, '2020-08-29': 'Too cloudy.', '2020-08-30': 3999.349569219829, '2020-08-31': 'Too cloudy.', '2020-09-01': 'Too cloudy.', '2020-09-02': 7082.840768113239, '2020-09-04': 'Too cloudy.', '2020-09-03': 'Too cloudy.', '2020-09-05': 'Too cloudy.', '2020-09-06': 3804.7631512784674, '2020-09-07': 4404.972897626298, '2020-09-09': 'Too cloudy.', '2020-09-08': 4792.396147471069, '2020-09-10': 4636.637770639921, '2020-09-11': 4108.221173168408, '2020-09-12': 5052.726613949424, '2020-09-13': 4870.492433329378, '2020-09-14': 'Too cloudy.', '2020-09-15': 'Too cloudy.', '2020-09-16': 'Too cloudy.', '2020-09-17': 'Too cloudy.', '2020-09-18': 2218.4663884015786, '2020-09-19': 'Too cloudy.', '2020-09-20': 2323.427012746887, '2020-09-21': 4598.338377157309, '2020-09-22': 3466.37126705203, '2020-09-23': 4578.797531328224, '2020-09-24': 4011.839579954459, '2020-09-25': 4654.024282718794, '2020-09-26': 'Too cloudy.', '2020-09-27': 7547.206540737594, '2020-09-28': 'Too cloudy.', '2020-09-29': 7006.630620502277, '2020-09-30': 4520.084491391584, '2020-10-01': 6013.600427716299, '2020-10-02': 4266.575547227999, '2020-10-03': 5121.750157264068, '2020-10-04': 4613.245211890603, '2020-10-05': 4156.7772146425705, '2020-10-06': 'Too cloudy.', '2020-10-07': 4853.266712559911, '2020-10-08': 7327.434646032423, '2020-10-09': 4586.626254892242, '2020-10-10': 2832.9385915554003, '2020-10-11': 4539.949078732267, '2020-10-12': 4342.9523313184645, '2020-10-13': 'Too cloudy.', '2020-10-14': 'Too cloudy.', '2020-10-15': 'Too cloudy.', '2020-10-16': 2690.976499049097, '2020-10-17': 4384.767994408241, '2020-10-18': 3042.093185515553, '2020-10-19': 'Too cloudy.', '2020-10-20': 'Too cloudy.', '2020-10-21': 'Too cloudy.', '2020-10-22': 'Too cloudy.', '2020-10-23': 4516.993091376657, '2020-10-24': 4573.561308243081, '2020-10-25': 4163.941665201333, '2020-10-26': 3661.0716346110926, '2020-10-27': 4568.28194232168, '2020-10-28': 'Too cloudy.', '2020-10-29': 'Too cloudy.', '2020-10-30': 4492.089781554009, '2020-10-31': 5009.1803327910475, '2020-11-01': 4494.526243731329, '2020-11-03': 'Too cloudy.', '2020-11-04': 4506.598966204573, '2020-11-05': 4190.50344487792, '2020-11-06': 'Too cloudy.', '2020-11-07': 4525.877496927874, '2020-11-08': 'Too cloudy.', '2020-11-09': 'Too cloudy.', '2020-11-10': 3414.912951677782, '2020-11-11': 4249.456207808294, '2020-11-12': 3733.945525504342, '2020-11-13': 3582.7312213418068, '2020-11-14': 'Too cloudy.', '2020-11-15': 4404.208231377214, '2020-11-16': 3866.544624680228, '2020-11-17': 3933.8607183735444, '2020-11-18': 'Too cloudy.', '2020-11-19': 'Too cloudy.', '2020-11-20': 3156.41657783258, '2020-11-21': 4476.519453005425, '2020-11-22': 2847.0320065817255, '2020-11-23': 3382.1639244933413, '2020-11-24': 'Too cloudy.', '2020-11-25': 4373.074730065389, '2020-11-26': 4063.59739541511, '2020-11-27': 4425.0963892276195, '2020-11-28': 'Too cloudy.', '2020-11-29': 4182.167156262863, '2020-11-30': 4149.612711602846, '2020-12-01': 4018.271422455172, '2020-12-02': 4397.067213836416, '2020-12-03': 4277.955387572472, '2020-12-05': 4165.154961940486, '2020-12-04': 4390.768405126966, '2020-12-06': 'Too cloudy.', '2020-12-07': 'Too cloudy.', '2020-12-08': 4223.75971999435, '2020-12-09': 4189.4206309373685, '2020-12-10': 'Too cloudy.', '2020-12-11': 4228.780296877149, '2020-12-12': 'Too cloudy.', '2020-12-13': 'Too cloudy.', '2020-12-14': 3946.0451302732686, '2020-12-15': 4041.003910418382, '2020-12-16': 'Too cloudy.', '2020-12-17': 'Too cloudy.', '2020-12-18': 4015.201585931147, '2020-12-19': 'Too cloudy.', '2020-12-20': 3853.8169277739207, '2020-12-21': 3785.828260062185, '2020-12-22': 3920.5632668648473, '2020-12-23': 3905.2186794192676, '2020-12-24': 3065.9374756458233, '2020-12-25': 3016.9341720601187, '2020-12-26': 'Too cloudy.', '2020-12-27': 'Too cloudy.', '2020-12-28': 'Too cloudy.', '2020-12-29': 2796.847223993264, '2020-12-30': 'Too cloudy.'}
data2021 = {'2021-01-15': 'Too cloudy.', '2021-01-16': 'Too cloudy.', '2021-01-17': 'Too cloudy.', '2021-01-19': 'Too cloudy.', '2021-01-18': 'Too cloudy.', '2021-01-21': 'Too cloudy.', '2021-01-20': 'Too cloudy.', '2021-01-22': 2654.255591201196, '2021-01-23': 2995.83321586727, '2021-01-24': 'Too cloudy.', '2021-01-25': 'Too cloudy.', '2021-01-26': 'Too cloudy.', '2021-01-29': 2795.0976091356524, '2021-01-27': 2503.449271483013, '2021-01-28': 2579.322759136135, '2021-01-30': 3334.4955240290487, '2021-01-31': 'Too cloudy.', '2021-02-01': 873.9065058925789, '2021-02-02': 4495.691030135837, '2021-02-03': 1286.8664958696443, '2021-02-04': 1200.779080680896, '2021-02-05': 5127.56808211595, '2021-02-06': 2539.5748862455157, '2021-02-07': 4211.049856784046, '2021-02-08': 1587.6750681920487, '2021-02-09': 2008.7114462518575, '2021-02-10': 709.5208436448654, '2021-02-11': 2778.097180730487, '2021-02-12': 929.6334142155852, '2021-02-13': 917.2478194889226, '2021-02-15': 'Too cloudy.', '2021-02-14': 4905.957333963093, '2021-02-16': 872.4660039846304, '2021-02-17': 1294.393133001706, '2021-02-18': 853.8543058318684, '2021-02-19': 730.6943805080497, '2021-02-20': 787.9342432782033, '2021-02-21': 686.3907584657086, '2021-02-22': 768.3008132692461, '2021-02-23': 'Too cloudy.', '2021-02-24': 'Too cloudy.', '2021-02-25': 1169.0980285083365, '2021-02-26': 2938.5016755201013, '2021-02-27': 'Too cloudy.', '2021-02-28': 'Too cloudy.', '2021-03-02': 769.878249173023, '2021-03-04': 634.4758556839004, '2021-03-01': 'Too cloudy.', '2021-03-03': 4995.828936458315, '2021-03-05': 2823.621085595843, '2021-03-06': 'Too cloudy.', '2021-03-07': 'Too cloudy.', '2021-03-08': 677.4793126819696, '2021-03-09': 3085.185042039111, '2021-03-10': 'Too cloudy.', '2021-03-11': 333.5760326707451, '2021-03-12': 1274.6004934604305, '2021-03-13': 'Too cloudy.', '2021-03-15': 5956.838909161205, '2021-03-14': 477.7192003508891, '2021-03-16': 432.7232352269526, '2021-03-17': 'Too cloudy.', '2021-03-18': 'Too cloudy.', '2021-03-19': 413.12761590003623, '2021-03-20': 2383.202348335021, '2021-03-21': 5825.30021128716, '2021-03-22': 1142.707964604985, '2021-03-23': 329.4879869686144, '2021-03-24': 'Too cloudy.', '2021-03-25': 'Too cloudy.', '2021-03-26': 'Too small.', '2021-03-27': 771.9602150022303, '2021-03-28': 466.28137930584455, '2021-03-29': 'Too cloudy.', '2021-03-30': 1562.7562974536502, '2021-03-31': 'Too cloudy.', '2021-04-01': 1897.7989345903156, '2021-04-03': 'Too cloudy.', '2021-04-02': 'Too cloudy.', '2021-04-04': 'Too cloudy.', '2021-04-05': 928.8693801444078, '2021-04-06': 'Too small.', '2021-04-07': 'Too small.', '2021-04-08': 'Too cloudy.', '2021-04-09': 'Too cloudy.', '2021-04-10': 630.76916430365, '2021-04-11': 'Too small.', '2021-04-12': 1535.5228972807897, '2021-04-13': 'Too cloudy.', '2021-04-14': 1414.258258652683, '2021-04-15': 'Too cloudy.', '2021-04-16': 'Too cloudy.', '2021-04-17': 'Too cloudy.', '2021-04-18': 'Too small.', '2021-04-19': 'Too cloudy.', '2021-04-20': 'Too small.', '2021-04-21': 'Too small.', '2021-01-14': 'Too cloudy.', '2021-01-13': 'Too cloudy.', '2021-01-12': 'Too cloudy.', '2021-01-11': 'Too cloudy.', '2021-01-10': 3209.9621622938766, '2021-01-09': 2788.205517285155, '2021-01-08': 'Too cloudy.', '2021-01-07': 2601.975341509975, '2021-01-06': 3431.5641380459174, '2021-01-05': 'Too cloudy.', '2021-01-04': 2820.427226747997, '2021-01-03': 'Too cloudy.', '2021-01-02': 'Too cloudy.', '2021-01-01': 2738.1934510831575}

#hand picked data
cleardates20 =['2020-01-15','2020-02-07','2020-02-25','2020-02-26','2020-03-22','2020-03-29','2020-04-09','2020-04-28','2020-05-28','2020-05-30','2020-06-15','2020-06-24','2020-07-14','2020-08-09','2020-08-24','2020-09-23','2020-10-11','2020-11-21','2020-11-30','2020-12-02','2020-12-11','2020-12-14','2020-12-24','2020-12-29','2020-12-31']
cleardata20 = [5605.867054748105,5718.880572004744,5356.204494016931,5282.54868945124,5248.619220428137,5325.483523561918,5377.414039670243,5340.442434208529,4965.275713804452,5079.820806020768,5123.097158208356,5074.976825609678,4866.498306857028,4634.825126872614,4883.642924131782,4578.797531594165,4539.949078997304,4476.519453261412,4149.612711840041,4397.067214084962,4228.780297102061,3946.045130475517,3065.9374757961423,2796.847224117944,2728.086308730213]
cleardates21 = ['2021-01-07','2021-01-09','2021-02-04','2021-02-10','2021-02-12','2021-02-16','2021-02-19','2021-03-02','2021-03-08','2021-03-16']
cleardata21 = [2601.975341635795,2788.2055174196694,1200.7790807340602,700.2703865042745,929.6334142560121,872.4660040218836,730.6943805376444,769.878249204907,677.4793127071263,432.72323524380147]


#filter data and merge
dates = []
data = []
for date in data2020:
    if type(data2020[date]) == float:
        if data2020[date]<=6000:
            if data2020[date]>=2000:
                dates.append(date)
                data.append(data2020[date])
    else:
        pass
dates.remove('2020-11-02')
data.remove(4051.092405713962)
dates21 = []
data21 = []
for date in data2021:
    if type(data2021[date]) == float:
        if data2021[date]<=4000:
            dates21.append(date)
            data21.append(data2021[date])
    else:
        pass
for date in dates21:
    dates.append(date)
for datum in data21:
    data.append(datum)

#create basic scatterplot
tnrfont = {'fontname':'Ubuntu'}

date = pd.to_datetime(dates)
DF = pd.DataFrame()
DF['value'] = data
DF = DF.set_index(date)
plt.scatter(DF.index, DF['value'], color='blue', s=16)
plt.gcf().autofmt_xdate()
plt.xticks(**tnrfont)
plt.yticks(**tnrfont)

plt.xlabel('Date',**tnrfont)
plt.ylabel('Area km2',**tnrfont)
plt.ylim(0,6000)
 
plt.title('A68a Area',**tnrfont, fontsize = 15)
 
plt.show()

#create comparison scatterplot
tnrfont = {'fontname':'Ubuntu'}

date = pd.to_datetime(dates)
DF = pd.DataFrame()
DF['value'] = data
DF = DF.set_index(date)
plt.scatter(DF.index, DF['value'], color='blue', s=16)
plt.gcf().autofmt_xdate()
plt.xticks(**tnrfont)
plt.yticks(**tnrfont)

date = pd.to_datetime(knowndates)
DF = pd.DataFrame()
DF['value'] = knowndata
DF = DF.set_index(date)
plt.scatter(DF.index, DF['value'], color='red', s=180, marker = 'x',linewidths=3)
plt.gcf().autofmt_xdate()

date = pd.to_datetime(cleardates20)
DF = pd.DataFrame()
DF['value'] = cleardata20
DF = DF.set_index(date)
plt.plot(DF,color='red',marker='o', markerfacecolor='red', markersize=2)
plt.gcf().autofmt_xdate()
date = pd.to_datetime(cleardates21)
DF = pd.DataFrame()
DF['value'] = cleardata21
DF = DF.set_index(date)
plt.plot(DF,color='red',marker='o', markerfacecolor='red', markersize=.5)
plt.gcf().autofmt_xdate()

plt.xlabel('Date',**tnrfont)
plt.ylabel('Area km2',**tnrfont)
plt.ylim(0,6000)
 
plt.title('A68a Area',**tnrfont, fontsize = 15)
 
plt.show()

