import argparse

args = argparse.ArgumentParser(description='manual to this script')
args.add_argument('--contigs', type=str, default='test_contigs.fa')
args.add_argument('--len', type=int, default=8000)
args.add_argument('--gpus', type=int, default=0)
args.add_argument('--model', type=str, default='pretrain', choices=['pretrain', 'retrain'], help='predicting host with pretrained parameters or retrained paramters (default pretrain)')
args.add_argument('--lr', default=0.001, type=float, help='Learning rate of training')
args.add_argument('--weight-decay', default=0, type=float, help='Weight Decay of optimizer')
args.add_argument('--topk',  type=int, default=1, help='topK')
args.add_argument('--t',  type=float, default=0.98, help='The confident threshold for predicting virus, the higier the threshold the higher the precision. (default 0.98)')
args.add_argument('--epochs', type=int, default=4000, help='Epochs to train')
args.add_argument('--batch-size', type=int, default=512, help='Batch size of Training')
args.add_argument('--GATf', type=int, default=512)
args.add_argument('--GATh', type=int, default=4)
args.add_argument('--head', type=int, default=4)


inputs = args.parse_args()
print(inputs)