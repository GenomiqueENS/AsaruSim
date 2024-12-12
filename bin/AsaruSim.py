import argparse
import version
from badread_caller import badread_caller, setup_badread_command
from template_maker import template_maker, setup_template_parameters
from PCR_amplificator import PCR_amplificator, setup_PCR_parameters
from truncation_estimator import truncation_estimator, setup_truncation_parameters

def setup_parent_parser():
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', action='store_true', help='Enable debug mode')
    return parent_parser


def main():
    parent_parser = setup_parent_parser()

    main_parser = argparse.ArgumentParser()
    main_parser.add_argument('--version', action='version', version=version.__version__)

    subparsers = main_parser.add_subparsers(title="command", dest="command", required=True)

    subparsers.add_parser('call_badread', parents=[setup_badread_command(parent_parser)], 
                                           add_help=False)
    subparsers.add_parser('template_maker', parents=[setup_template_parameters(parent_parser)], 
                                             add_help=False)
    subparsers.add_parser('PCR', parents=[setup_PCR_parameters(parent_parser)], 
                                             add_help=False)
    subparsers.add_parser('truncation_estimator', parents=[setup_truncation_parameters(parent_parser)], 
                                             add_help=False)

    args = main_parser.parse_args()

    if args.command == 'call_badread':
        badread_caller(args)
    if args.command == 'template_maker':
        template_maker(args)
    if args.command == 'PCR':
        PCR_amplificator(args)
    if args.command == 'truncation_estimator':
        truncation_estimator(args)
    

if __name__ == "__main__":
    main()
