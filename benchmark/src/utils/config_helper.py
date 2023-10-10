import yaml

def get_data_config(dataset_name = None, overrides = None):
  with open('data/datasets.yaml', 'r') as file:
    config = yaml.safe_load(file)
    if dataset_name is None and overrides is None:
      pass
    else:
      if dataset_name not in config.keys():
        raise ValueError(f"No data config defined for {dataset_name}.\nThe options are {', '.join(config.keys())}")
      config = config[dataset_name]
    
    if overrides is None or overrides == "":
      return config
    else:
      if type(overrides) in [list, tuple]:
        return update_config(config, overrides)
      else:
        return update_config(config, [overrides])
        
  
def update_config(config, update):
    lhs_rhs = [item.split("=", 1) for item in update]
    for pair in lhs_rhs:
        lhs_elem = pair[0].split(".")
        nested_dict = config
        for elem in lhs_elem[:-1]:
            nested_dict = nested_dict.get(elem, {})
        nested_dict[lhs_elem[-1]] = yaml.safe_load(pair[1])
    return config
