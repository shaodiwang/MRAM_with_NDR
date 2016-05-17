sed s/num_trials/65536/g start_simulation.template > start_simulation.py
python start_simulation.py
sed s/num_trials/1048576/g start_simulation.template > start_simulation.py
python start_simulation.py
sed s/num_trials/10485760/g start_simulation.template > start_simulation.py
python start_simulation.py
sed s/num_trials/104857600/g start_simulation.template > start_simulation.py
python start_simulation.py
sed s/num_trials/419430400/g start_simulation.template > start_simulation.py
python start_simulation.py
