"""Add gate table

Revision ID: 06df70ede08a
Revises: fa29b07fd05f
Create Date: 2020-06-10 17:37:17.185233

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB, ARRAY


# revision identifiers, used by Alembic.
revision = '06df70ede08a'
down_revision = 'fa29b07fd05f'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'gate',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('dataset_id', sa.Integer(), sa.ForeignKey("dataset.id", ondelete="CASCADE"), nullable=False, index=True),
        sa.Column('name', sa.String()),
        sa.Column('description', sa.String()),
        sa.Column('acquisition_ids', ARRAY(sa.Integer())),
        sa.Column('indices', ARRAY(sa.Integer())),
        sa.Column('cell_ids', ARRAY(sa.Integer())),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('gate')
