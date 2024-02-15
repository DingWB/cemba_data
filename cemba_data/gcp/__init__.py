from .demultiplex import *
try:
	import fire
except:
	pip_path = os.path.join(os.path.dirname(sys.executable), 'pip')
	os.system(f"{pip_path} install fire")
	import fire

def main():
	fire.core.Display = lambda lines, out: print(*lines, file=out)
	# fire.Fire()
	fire.Fire({
		"prepare_demultiplex":prepare_demultiplex,
		"get_demultiplex_skypilot_yaml":get_demultiplex_skypilot_yaml,
		'run_demultiplex':run_demultiplex,
		'prepare_mapping':prepare_mapping,
		'run_mapping':run_mapping,
		'yap_pipeline':yap_pipeline,
		'check_demultiplex':check_demultiplex,
		'cell_qc':cell_qc,
	})

if __name__=="_main__":
	main()