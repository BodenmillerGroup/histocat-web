"""Add share table

Revision ID: 35b87e175854
Revises: 9d5cad983bce
Create Date: 2019-07-27 21:56:03.229341

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import ARRAY


# revision identifiers, used by Alembic.
revision = '35b87e175854'
down_revision = '9d5cad983bce'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'share',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('user_id', sa.Integer(), sa.ForeignKey("user.id", ondelete="CASCADE"), nullable=False, index=True),
        sa.Column('experiment_id', sa.Integer(), sa.ForeignKey("experiment.id", ondelete="CASCADE"), nullable=False, index=True),
        sa.Column('permissions', ARRAY(sa.String(64))),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )
    op.create_unique_constraint(
        'uq_share_user_id_experiment_id',
        'share',
        ['user_id', 'experiment_id']
    )


def downgrade():
    op.drop_table('share')
    op.drop_unique_constraint('uq_share_user_id_experiment_id')
