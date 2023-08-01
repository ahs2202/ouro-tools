import numpy as np


# In[ ]:


dict_a2b = dict()  # an empty dictionary for mapping


# In[ ]:


def Remove_version_info(ID):
    return ID.split(".")[0]


# In[ ]:


def Retrive_Length(entry):
    return len(entry)


class Map(object):
    def __init__(self, dict_a2b):
        self.dict_a2b = dict_a2b

    def a2b(self, a):
        if a in self.dict_a2b:
            return self.dict_a2b[a]
        else:
            return np.nan

    def a2b_if_mapping_available_else_Map_a2a(self, a):
        if a in self.dict_a2b:
            return self.dict_a2b[a]
        else:
            return a
