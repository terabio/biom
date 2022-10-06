import torch
from torch import nn


class ResidualBlock(nn.Module):
    def __init__(self, in_channels, out_channels, kernel_size, padding):
        super().__init__()
        self.seq = nn.Sequential(
            nn.Conv1d(in_channels, out_channels, kernel_size=kernel_size, padding=padding),
            nn.BatchNorm1d(out_channels),
            nn.ReLU()
        )
        self.adjust = nn.Sequential(
            nn.Conv1d(in_channels, out_channels, kernel_size=1),
            nn.BatchNorm1d(out_channels)
        ) if in_channels != out_channels else lambda x: x
        self.bn = nn.BatchNorm1d(out_channels)

    def forward(self, x):
        x = self.seq(x) + self.adjust(x)
        x = self.bn(x).relu()
        return x


class SEBlock(nn.Module):
    def __init__(self, inchannels):
        super().__init__()
        self.pooling = nn.AdaptiveAvgPool1d(1)
        self.fc = nn.Linear(inchannels, inchannels)

    def forward(self, x):
        chns = self.pooling(x).squeeze(-1)
        chns = self.fc(chns).unsqueeze(-1)
        return x + chns


class UNet(nn.Module):
    def __init__(self, filters_scale=1, poolsize=4):
        super().__init__()
        fs = filters_scale
        self.encoder_block1 = nn.Sequential(
            ResidualBlock(10, int(10 * fs), kernel_size=11, padding=5),
            ResidualBlock(int(10 * fs), int(12 * fs), kernel_size=11, padding=5),
            ResidualBlock(int(12 * fs), int(16 * fs), kernel_size=11, padding=5),
            # SEBlock(int(16 * mul)),
            # ResidualBlock(int(16 * mul), int(16 * mul), kernel_size=11, padding=5),
            ResidualBlock(int(16 * fs), int(16 * fs), kernel_size=11, padding=5)
        )
        self.pool1 = nn.MaxPool1d(kernel_size=poolsize)
        self.encoder_block2 = nn.Sequential(
            ResidualBlock(int(16 * fs), int(16 * fs), kernel_size=11, padding=5),
            ResidualBlock(int(16 * fs), int(32 * fs), kernel_size=11, padding=5),
            # SEBlock(int(32 * mul)),
            ResidualBlock(int(32 * fs), int(32 * fs), kernel_size=11, padding=5),
            ResidualBlock(int(32 * fs), int(32 * fs), kernel_size=11, padding=5),
        )
        self.pool2 = nn.MaxPool1d(kernel_size=poolsize)
        self.encoder_block3 = nn.Sequential(
            ResidualBlock(int(32 * fs), int(32 * fs), kernel_size=5, padding=2),
            ResidualBlock(int(32 * fs), int(64 * fs), kernel_size=5, padding=2),
            # SEBlock(int(64 * mul)),
            ResidualBlock(int(64 * fs), int(64 * fs), kernel_size=5, padding=2),
            ResidualBlock(int(64 * fs), int(64 * fs), kernel_size=5, padding=2),
        )
        self.pool3 = nn.MaxPool1d(kernel_size=poolsize)
        self.encoder_block4 = nn.Sequential(
            ResidualBlock(int(64 * fs), int(64 * fs), kernel_size=5, padding=2),
            ResidualBlock(int(64 * fs), int(128 * fs), kernel_size=5, padding=2),
            # SEBlock(int(128 * mul)),
            ResidualBlock(int(128 * fs), int(128 * fs), kernel_size=5, padding=2),
            ResidualBlock(int(128 * fs), int(128 * fs), kernel_size=5, padding=2)
        )
        self.upsample1 = nn.Upsample(scale_factor=poolsize)
        self.decoder_block1 = nn.Sequential(
            nn.Conv1d(int(128 * fs) + int(64 * fs), int(128 * fs), kernel_size=1),
            nn.BatchNorm1d(int(128 * fs)),
            nn.ReLU(),

            ResidualBlock(int(128 * fs), int(64 * fs), kernel_size=5, padding=2),
            # SEBlock(int(64 * mul)),
            ResidualBlock(int(64 * fs), int(64 * fs), kernel_size=5, padding=2),
            ResidualBlock(int(64 * fs), int(64 * fs), kernel_size=5, padding=2)
        )
        self.upsample2 = nn.Upsample(scale_factor=poolsize)
        self.decoder_block2 = nn.Sequential(
            nn.Conv1d(int(64 * fs) + int(32 * fs), int(64 * fs), kernel_size=1),
            nn.BatchNorm1d(int(64 * fs)),
            nn.ReLU(),

            ResidualBlock(int(64 * fs), int(32 * fs), kernel_size=5, padding=2),
            # SEBlock(int(32 * mul)),
            ResidualBlock(int(32 * fs), int(32 * fs), kernel_size=5, padding=2),
            ResidualBlock(int(32 * fs), int(32 * fs), kernel_size=5, padding=2)
        )
        self.upsample3 = nn.Upsample(scale_factor=poolsize)
        self.decoder_block3 = nn.Sequential(
            nn.Conv1d(int(32 * fs) + int(16 * fs), int(32 * fs), kernel_size=1),
            nn.BatchNorm1d(int(32 * fs)),
            nn.ReLU(),

            ResidualBlock(int(32 * fs), int(16 * fs), kernel_size=11, padding=5),
            # SEBlock(int(16 * mul)),
            ResidualBlock(int(16 * fs), int(16 * fs), kernel_size=11, padding=5),
            ResidualBlock(int(16 * fs), int(8 * fs), kernel_size=11, padding=5),
            ResidualBlock(int(8 * fs), int(8 * fs), kernel_size=11, padding=5),
            nn.Conv1d(int(8 * fs), 4, kernel_size=11, padding=5)
        )

    def forward(self, x):
        assert x.shape[1] == 10
        # down pass
        stage1 = self.encoder_block1(x)
        stage2 = self.encoder_block2(self.pool1(stage1))
        stage3 = self.encoder_block3(self.pool2(stage2))
        predict = self.encoder_block4(self.pool3(stage3))

        # up pass
        predict = self.upsample1(predict)
        predict = self.decoder_block1(torch.cat([predict, stage3], dim=1))
        predict = self.upsample2(predict)
        predict = self.decoder_block2(torch.cat([predict, stage2], dim=1))
        predict = self.upsample3(predict)
        predict = self.decoder_block3(torch.cat([predict, stage1], dim=1))

        N, _, L = x.shape
        assert (N, 4, L) == tuple(predict.shape), predict.shape
        return predict.sigmoid()
