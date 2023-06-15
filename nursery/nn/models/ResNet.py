from typing import List, Optional, Tuple

import torch.nn as nn
from torch import Tensor


def conv7(in_planes: int, out_planes: int, stride: int = 1) -> nn.Conv1d:
    return nn.Conv1d(in_planes, out_planes, kernel_size=7, stride=stride, bias=False, padding=3)


def conv1(in_planes: int, out_planes: int, stride: int = 1) -> nn.Conv1d:
    return nn.Conv1d(in_planes, out_planes, kernel_size=1, stride=stride, bias=False)


class ResNetBlock(nn.Module):
    def __init__(
            self,
            inplanes: int,
            planes: int,
            stride: int = 1,
            downsample: Optional[nn.Module] = None,
    ) -> None:
        super(ResNetBlock, self).__init__()
        # Both self.conv1 and self.downsample layers downsample the input when stride != 1
        self.conv1 = conv7(inplanes, planes, stride)
        self.bn1 = nn.BatchNorm1d(planes)
        self.relu = nn.ReLU(inplace=True)
        self.conv2 = conv7(planes, planes)
        self.bn2 = nn.BatchNorm1d(planes)
        self.downsample = downsample
        self.stride = stride

    def forward(self, x: Tensor) -> Tensor:
        identity = x

        out = self.conv1(x)
        out = self.bn1(out)
        out = self.relu(out)

        out = self.conv2(out)
        out = self.bn2(out)

        if self.downsample is not None:
            identity = self.downsample(x)

        out += identity
        out = self.relu(out)

        return out


class ResNet(nn.Module):
    def __init__(
            self,
            inplanes: int,
            planes: List[int],
            strides: List[int],
            classes: int,
            dropout: float = 0.5,
    ) -> None:
        super(ResNet, self).__init__()
        self.backbone, depth = self._make_backbone(inplanes, planes, strides)
        self.dropout = nn.Dropout(p=dropout)
        self.fc = nn.Linear(depth, classes)

    def _make_backbone(self, inplanes: int, planes: List[int], strides: List[int]) -> Tuple[nn.Sequential, int]:
        layers = []
        depth = inplanes
        for pln, strd in zip(planes, strides):
            downsample = nn.Sequential(
                conv1(depth, pln, strd),
                nn.BatchNorm1d(pln)
            )
            layers.append(ResNetBlock(depth, pln, strd, downsample))
            depth = pln

        layers.append(nn.AdaptiveMaxPool1d(1))
        layers.append(nn.Flatten())
        layers = nn.Sequential(*layers)
        # # Fancy initialization to improve performance
        # for m in layers.modules():
        #     raise Exception
        #     # Must be done inside ResNetBlock
        #     if isinstance(m, nn.Conv2d):
        #         nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
        #     elif isinstance(m, (nn.BatchNorm2d, nn.GroupNorm)):
        #         nn.init.constant_(m.weight, 1)
        #         nn.init.constant_(m.bias, 0)
        #     elif isinstance(m, resblock.ResNetBlock):
        #         nn.init.constant_(m.bn2.weight, 0)
        return layers, depth

    def forward(self, indata: Tensor) -> Tensor:
        x = self.backbone(indata)
        x = self.dropout(x)
        x = self.fc(x)
        return x
