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
		'gcp_yap_pipeline':gcp_yap_pipeline,
	})

if __name__=="_main__":
	main()