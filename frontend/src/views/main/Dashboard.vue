<template>
  <v-flex>
    <v-combobox
      v-model="model"
      :filter="filter"
      :hide-no-data="!search"
      :items="items"
      :search-input.sync="search"
      hide-selected
      label="Search for an option"
      multiple
      small-chips
      solo
      class="mt-3 mx-3"
    >
      <template v-slot:no-data>
        <v-list-tile>
          <span class="subheading">Create</span>
          <v-chip
            :color="`${colors[nonce - 1]} lighten-3`"
            label
            small
          >
            {{ search }}
          </v-chip>
        </v-list-tile>
      </template>
      <template v-slot:selection="{ item, parent, selected }">
        <v-chip
          v-if="item === Object(item)"
          :color="`${item.color} lighten-3`"
          :selected="selected"
          label
          small
        >
        <span class="pr-2">
          {{ item.text }}
        </span>
          <v-icon
            small
            @click="parent.selectItem(item)"
          >mdi-close
          </v-icon>
        </v-chip>
      </template>
      <template v-slot:item="{ index, item }">
        <v-list-tile-content>
          <v-text-field
            v-if="editing === item"
            v-model="editing.text"
            autofocus
            flat
            background-color="transparent"
            hide-details
            solo
            @keyup.enter="edit(index, item)"
          ></v-text-field>
          <v-chip
            v-else
            :color="`${item.color} lighten-3`"
            dark
            label
            small
          >
            {{ item.text }}
          </v-chip>
        </v-list-tile-content>
        <v-spacer></v-spacer>
        <v-list-tile-action @click.stop>
          <v-btn
            icon
            @click.stop.prevent="edit(index, item)"
          >
            <v-icon>{{ editing !== item ? 'mdi-pencil' : 'mdi-check' }}</v-icon>
          </v-btn>
        </v-list-tile-action>
      </template>
    </v-combobox>

    <masonry
      :cols="3"
      :gutter="30"
    >
      <ExperimentCard
        v-for="experiment in experiments"
        :experiment="experiment"
        :key="experiment.id"
      />
    </masonry>
  </v-flex>
</template>

<script lang="ts">
  import { Component, Vue, Watch } from 'vue-property-decorator';
  import { readAdminExperiments } from '@/modules/experiment/getters';
  import { dispatchGetExperiments } from '@/modules/experiment/actions';
  import ExperimentCard from '@/components/ExperimentCard.vue';

  @Component({
    components: { ExperimentCard },
  })
  export default class Dashboard extends Vue {
    activator = null;
    attach = null;
    colors = ['green', 'purple', 'indigo', 'cyan', 'teal', 'orange'];
    editing = null;
    index = -1;
    items = [
      { header: 'Select an option or create one' },
      {
        text: 'Foo',
        color: 'blue',
      },
      {
        text: 'Bar',
        color: 'red',
      },
    ];
    nonce = 1;
    menu = false;
    model = [
      {
        text: 'Foo',
        color: 'blue',
      },
    ];
    x = 0;
    search = null;
    y = 0;

    get experiments() {
      return readAdminExperiments(this.$store);
    }

    @Watch('model')
    onModelChange(
      val, prev,
    ) {
      if (val.length === prev.length) {
        return;
      }

      this.model = val.map(v => {
        if (typeof v === 'string') {
          v = {
            text: v,
            color: this.colors[this.nonce - 1],
          };

          this.items.push(v);

          this.nonce++;
        }

        return v;
      });
    }

    edit(index: number, item) {
      if (!this.editing) {
        this.editing = item;
        this.index = index;
      } else {
        this.editing = null;
        this.index = -1;
      }
    }

    filter(item, queryText: string, itemText: string) {
      if (item.header) {
        return false;
      }

      const hasValue = val => val != null ? val : '';

      const text = hasValue(itemText);
      const query = hasValue(queryText);

      return text.toString()
        .toLowerCase()
        .indexOf(query.toString().toLowerCase()) > -1;
    }

    async mounted() {
      await dispatchGetExperiments(this.$store);
    }
  }
</script>
