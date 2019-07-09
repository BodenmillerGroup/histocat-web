"""Add aquisition table

Revision ID: b182e0878258
Revises: c537dg62e8je
Create Date: 2019-04-10 15:50:16.499555

"""
import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import JSONB

# revision identifiers, used by Alembic.
revision = 'b182e0878258'
down_revision = 'c537dg62e8je'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'acquisition',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('roi_id', sa.Integer(), sa.ForeignKey("roi.id", ondelete="CASCADE"), index=True),
        sa.Column('description', sa.String()),
        sa.Column('metaname', sa.String(4096)),
        sa.Column('original_id', sa.Integer(), index=True),
        sa.Column('order_number', sa.Integer()),
        sa.Column('ablation_power', sa.Float()),
        sa.Column('ablation_distance_between_shots_x', sa.Float()),
        sa.Column('ablation_distance_between_shots_y', sa.Float()),
        sa.Column('ablation_frequency', sa.Float()),
        sa.Column('signal_type', sa.String(64)),
        sa.Column('dual_count_start', sa.Integer()),
        sa.Column('data_start_offset', sa.BigInteger()),
        sa.Column('data_end_offset', sa.BigInteger()),
        sa.Column('start_timestamp', sa.DateTime()),
        sa.Column('end_timestamp', sa.DateTime()),
        sa.Column('after_ablation_image_end_offset', sa.BigInteger()),
        sa.Column('after_ablation_image_start_offset', sa.BigInteger()),
        sa.Column('before_ablation_image_end_offset', sa.BigInteger()),
        sa.Column('before_ablation_image_start_offset', sa.BigInteger()),
        sa.Column('roi_start_x_pos_um', sa.Float()),
        sa.Column('roi_start_y_pos_um', sa.Float()),
        sa.Column('roi_end_x_pos_um', sa.Float()),
        sa.Column('roi_end_y_pos_um', sa.Float()),
        sa.Column('movement_type', sa.String(64)),
        sa.Column('segment_data_format', sa.String(64)),
        sa.Column('value_bytes', sa.Integer()),
        sa.Column('max_y', sa.Integer()),
        sa.Column('max_x', sa.Integer()),
        sa.Column('plume_start', sa.Integer()),
        sa.Column('plume_end', sa.Integer()),
        sa.Column('template', sa.String()),
        sa.Column('meta', JSONB()),
        sa.Column('location', sa.String(4096)),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('acquisition')
