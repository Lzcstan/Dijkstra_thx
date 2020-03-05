#对于word的定义：采取OfficeWord中的标准，每一连续字母串、数字都算作一个word，对于#等排版符号，在读取时预先删除；文件名称和网址每一个算作一个word
def word_count(fn):
    import re
    global words
    words = {}
    r = re.compile(r"[,!*.#]")
    with open(fn,"r") as f:
        for line in f:
            for word in r.sub("",line.strip()).split(" "):
             if word in words:
                    words[word] += 1
             words.setdefault(word,1)
    del words['']
    print('The dictionary of readme.md is',words,'\n',sum(words.values()))
word_count('C:/Users/lzcstan/Desktop/git_test/SES2020spring/unit2/readme.md')
