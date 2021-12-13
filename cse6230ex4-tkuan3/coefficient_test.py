import json

with open("coefficient_test.json", "r") as read_file:
  data = json.load(read_file)
  v = data["estimated diffusivities"][1]
  if 0.799815 < v and v < 0.876264:
    print("coefficient test passed")
  else:
    print(f"coefficient test not passed: coefficient {v} is not in the acceptable range (0.799815, 0.876264)")
    raise RuntimeError


