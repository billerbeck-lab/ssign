# Vendored from https://github.com/zhengdd0422/PLM-Effector (CC-BY 3.0).
# See the package __init__.py for attribution. Upstream imports of
# torch_geometric, `from utils import pool_features`, and the unused
# transformers symbols have been pruned — none of the classes defined
# here reference those names at runtime.
import math

import torch
import torch.nn as nn
import torch.nn.functional as F


class SimpleMLP(nn.Module):
    def __init__(self, input_dim, hidden_dim, hidden_layer=2, dropout=0.5):
        super(SimpleMLP, self).__init__()
        self.dp = hidden_layer
        self.dropout = nn.Dropout(dropout)
        self.fc1 = nn.Linear(input_dim, hidden_dim)
        self.fc2 = nn.Linear(hidden_dim, hidden_dim)
        self.fc_out = nn.Linear(hidden_dim, 1)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        x = torch.relu(self.fc1(x))
        x = self.dropout(x)

        if self.dp == 2:
            x = torch.relu(self.fc2(x))
            x = self.dropout(x)

        x = self.fc_out(x)
        # x = self.sigmoid(x)

        return x


class SimpleCNN1D(nn.Module):
    def __init__(
        self,
        input_dim,
        cnn_depth,
        num_filters=64,
        kernel_size=3,
        maxpool=5,
        fc_hidden_dim=128,
        num_classes=1,
        dropout=0.5,
    ):
        super(SimpleCNN1D, self).__init__()
        self.cnn_dp = cnn_depth
        self.conv1 = nn.Conv1d(
            in_channels=input_dim,
            out_channels=num_filters,
            kernel_size=kernel_size,
            padding=1,
        )
        self.conv2 = nn.Conv1d(
            in_channels=num_filters,
            out_channels=num_filters,
            kernel_size=kernel_size,
            padding=1,
        )
        self.pool = nn.MaxPool1d(kernel_size=maxpool)
        self.fc = nn.Linear(num_filters, fc_hidden_dim)  # 直接定义好，不放在 forward
        self.fc_out = nn.Linear(fc_hidden_dim, num_classes)  # 输出层
        self.dropout = nn.Dropout(dropout)

    def forward(self, x):
        # 转换输入形状：变成 (batch_size, feature_dim, seq_len)
        x = x.permute(0, 2, 1)  # 重要！Conv1d 期望的是 (batch, in_channels, seq_len)

        x = torch.relu(self.conv1(x))
        x = self.pool(x)  # 池化

        if self.cnn_dp == 2:
            x = torch.relu(self.conv2(x))
            x = self.pool(x)  # 再次池化

        x = torch.mean(x, dim=-1)  # 全局平均池化
        x = self.fc(x)
        x = torch.relu(x)
        x = self.dropout(x)

        x = self.fc_out(x)

        return x


class CNNWithAttention(nn.Module):
    def __init__(
        self,
        input_dim,
        cnn_depth,
        num_filters=64,
        kernel_size=3,
        maxpool=5,
        fc_hidden_dim=128,
        num_classes=1,
        dropout=0.5,
    ):
        super(CNNWithAttention, self).__init__()

        # 1. 保留原有CNN部分
        self.cnn_dp = cnn_depth
        self.conv1 = nn.Conv1d(
            in_channels=input_dim,
            out_channels=num_filters,
            kernel_size=kernel_size,
            padding=1,
        )
        self.conv2 = nn.Conv1d(
            in_channels=num_filters,
            out_channels=num_filters,
            kernel_size=kernel_size,
            padding=1,
        )
        self.pool = nn.MaxPool1d(kernel_size=maxpool)

        # 2. 新增Attention模块
        self.attention = nn.MultiheadAttention(
            embed_dim=num_filters,  # 输入维度与CNN输出通道一致
            num_heads=4,  # 4头注意力
            dropout=dropout,
            batch_first=True,  # 输入格式为(batch, seq, features)
        )
        self.attn_layer_norm = nn.LayerNorm(num_filters)

        # 3. 分类头（保持不变）
        self.fc = nn.Linear(num_filters, fc_hidden_dim)
        self.fc_out = nn.Linear(fc_hidden_dim, num_classes)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x):
        # 输入形状: (batch_size, seq_len, input_dim)
        x = x.permute(0, 2, 1)  # -> (batch, input_dim, seq_len)

        # 原始CNN部分
        x = F.relu(self.conv1(x))
        x = self.pool(x)
        if self.cnn_dp == 2:
            x = F.relu(self.conv2(x))
            x = self.pool(x)

        # CNN输出形状: (batch, num_filters, seq_len)
        x = x.permute(0, 2, 1)  # -> (batch, seq_len, num_filters) 适应Attention

        # Attention部分
        attn_out, _ = self.attention(x, x, x)  # 自注意力
        x = self.attn_layer_norm(x + attn_out)  # 残差连接

        # 全局平均池化 + 分类头
        x = x.mean(dim=1)  # (batch, num_filters)
        x = F.relu(self.fc(x))
        x = self.dropout(x)
        return self.fc_out(x)  # 输出logits


class PositionalEncoding(nn.Module):
    """位置编码（可选）"""

    def __init__(self, d_model: int, max_len: int = 1024):
        super().__init__()
        position = torch.arange(max_len).unsqueeze(1)
        div_term = torch.exp(
            torch.arange(0, d_model, 2) * (-math.log(10000.0) / d_model)
        )
        pe = torch.zeros(max_len, d_model)
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        self.register_buffer("pe", pe)  # (max_len, d_model)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """在输入特征中加入位置编码"""
        return x + self.pe[: x.size(1)]  # 自动适配实际序列长度


class AttentionOnlyModel(nn.Module):
    def __init__(
        self,
        input_dim,
        max_len,
        hidden_dim,
        num_heads=8,
        dropout=0.1,
        num_classes=1,
        model_type="esm1",
    ):
        """
        Args:
            input_dim: 输入的特征维度
            max_len: 序列的最大长度
            num_heads: 注意力头的数量
            dropout: Dropout的比例
            num_classes: 分类类别数，二分类为1
            model_type: 指定使用的模型类型
        """
        super(AttentionOnlyModel, self).__init__()
        self.model_type = model_type

        # 位置编码（针对长序列）
        self.pos_encoder = PositionalEncoding(input_dim, max_len=max_len)

        # 注意力层
        self.attn = nn.MultiheadAttention(
            embed_dim=input_dim, num_heads=num_heads, batch_first=True
        )

        # 归一化层
        self.norm = nn.LayerNorm(input_dim)

        # Dropout层
        self.dropout = nn.Dropout(dropout)

        self.fc1 = nn.Linear(input_dim, hidden_dim)
        self.relu = nn.ReLU()
        # 输出层（根据二分类或多分类输出结果）
        self.fc_out = nn.Linear(hidden_dim, num_classes)

        # 需要遮蔽的标记
        self.cls_token_pos = 0  # ESM2的[CLS]位置
        self.eos_token_pos = -1  # ESM2的[EOS]位置

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: 输入特征张量，形状为(batch, seq_len, input_dim)
        """
        # 创建初始mask
        mask = torch.ones(x.size(1), dtype=torch.bool, device=x.device)

        # 根据模型类型来修改mask
        if self.model_type == "protT5":
            # 对于ProtT5模型，只遮蔽[EOS]（末尾token）
            mask[x.size(1) - 1] = False
        else:
            # 对于ESM等模型，遮蔽[CLS]和[EOS]
            mask[self.cls_token_pos] = False
            mask[self.eos_token_pos] = False

        # 将mask扩展为batch大小
        attn_mask = mask.unsqueeze(0).expand(x.size(0), -1)  # (batch, seq_len)

        # 对输入添加位置编码
        x = self.pos_encoder(x)

        # 注意力计算
        attn_output, _ = self.attn(
            x, x, x, key_padding_mask=~attn_mask if self.training else None
        )

        # 残差连接 + 层归一化
        x = self.norm(x + attn_output)

        # 聚合特征：全局平均池化
        x_pooled = x.mean(dim=1)

        # Dropout和输出
        x = self.dropout(x_pooled)
        x = self.fc1(x)
        x = self.relu(x)
        x = self.dropout(x)
        x = self.fc_out(x)

        return x
