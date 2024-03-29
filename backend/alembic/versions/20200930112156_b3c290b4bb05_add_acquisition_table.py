"""Add acquisition table

Revision ID: b3c290b4bb05
Revises: 69e599da5845
Create Date: 2020-09-30 11:21:56.231395

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB


# revision identifiers, used by Alembic.
revision = 'b3c290b4bb05'
down_revision = '69e599da5845'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'acquisition',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('slide_id', sa.Integer(), sa.ForeignKey("slide.id", ondelete="CASCADE"), index=True, nullable=False),
        sa.Column('origin_id', sa.Integer(), index=True),
        sa.Column('description', sa.String()),
        sa.Column('max_x', sa.Integer()),
        sa.Column('max_y', sa.Integer()),
        sa.Column('signal_type', sa.String()),
        sa.Column('segment_data_format', sa.String()),
        sa.Column('ablation_frequency', sa.Float()),
        sa.Column('ablation_power', sa.Float()),
        sa.Column('start_timestamp', sa.DateTime(timezone=True)),
        sa.Column('end_timestamp', sa.DateTime(timezone=True)),
        sa.Column('movement_type', sa.String()),
        sa.Column('ablation_distance_between_shots_x', sa.Float()),
        sa.Column('ablation_distance_between_shots_y', sa.Float()),
        sa.Column('template', sa.String()),
        sa.Column('roi_start_x_pos_um', sa.Float()),
        sa.Column('roi_start_y_pos_um', sa.Float()),
        sa.Column('roi_end_x_pos_um', sa.Float()),
        sa.Column('roi_end_y_pos_um', sa.Float()),
        sa.Column('has_before_ablation_image', sa.Boolean(), nullable=False, default=False),
        sa.Column('has_after_ablation_image', sa.Boolean(), nullable=False, default=False),
        sa.Column('is_valid', sa.Boolean(), nullable=False, default=True),
        sa.Column('channels', JSONB()),
        sa.Column('meta', JSONB()),
        sa.Column('location', sa.String(4096)),
    )


def downgrade():
    op.drop_table('acquisition')
